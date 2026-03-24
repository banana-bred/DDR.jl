import Plots

export plot_eigenph, plot_eigenph!
export plot_resonance, plot_resonance!
export overlay_resonances!

#################
#### HELPERS ####
#################
function _default_ylabel(ytype :: Symbol)
  ytype === :phase    && return "δ"
  ytype === :deriv    && return "dδ/dE"
  ytype === :absderiv && return "|dδ/dE|"
  error("Symbol ytype must be :phase, :deriv, or :absderiv")
end

function _default_record_label(record :: EigenphRecord)
  g = geom_index(record.path)
  return g === nothing  ? splitpath(record.path)[end] : "geom = $(g)"
end

################
#### PUBLIC ####
################

function plot_eigenph!(
    plt
  , data :: EigenphData
  ; cols=:all
  , selector=nothing
  , ytype :: Symbol = :phase
  , smoothwin :: Int = 1
  , label = :name
  , unitsE :: Symbol = :eV
  , kwargs...
  )
  colidx = resolve_cols(data; cols=cols, selector=selector)
  isempty(colidx) && error("No columns selected")
  Evals = unitsE === :au ? data.E :
          unitsE === :eV ? _convertE.(data.E, :au, :eV) :
          error("`unitsE` must be either :au or :eV")
  for c in colidx
    ys = _yseries(data, c; ytype=ytype, smoothwin=smoothwin)
    lbl = label isa Function ? label(data, c) :
      label === :name ? data.names[c] :
        label === false || label === :none || label === nothing ? false :
        string(label)
    (ytype === :deriv || ytype === :absderiv) && (ys = _convertInvE.(ys, :au, unitsE))
    Plots.plot!(plt, Evals, ys; label=lbl, kwargs...)
  end
  return plt
end

function plot_eigenph(
    data :: EigenphData
  ; cols=:all
  , selector=nothing
  , ytype :: Symbol = :phase
  , smoothwin :: Int = 1
  , label=:name
  , unitsE :: Symbol = :eV

  # -- general plot stuff
  , ylabel = nothing
  , xlabel = "E"

  # -- resonance overlay control
  , resonances=nothing  # nothing | :auto | Vector{<:Resonance}
  , res_detect :: NamedTuple = (;) # kwargs forwarded to find_resonances when resonances=:auto

  # -- overlay styling
  , res_color=:black
  , res_alpha :: Real=0.8
  , res_linewidth=0.8
  , res_linestyle=:solid
  , res_show_center :: Bool=true
  , res_markersize :: Real=3

  # -- other plotting kwargs
  , kwargs...
  )

  colidx = resolve_cols(data; cols=cols, selector=selector)
  isempty(colidx) && error("No columns selected !")

  ylbl = ylabel === nothing ? _default_ylabel(ytype) : ylabel
  plt = Plots.plot(; ylabel=ylbl, xlabel=xlabel)

  plot_eigenph!( plt
               , data
               ; cols=cols
               , selector=nothing
               , ytype=ytype
               , smoothwin=smoothwin
               , label=label
               , unitsE = unitsE
               , kwargs...
  )

  reslist = resolve_resonances(data, resonances, colidx, res_detect)
  reslist === nothing && return plt

  overlay_resonances!(plt, data, reslist
                    ; ytype=ytype
                    , smoothwin=smoothwin
                    , unitsE = unitsE
                    , cols_selected=colidx
                    , res_color=res_color
                    , res_alpha=res_alpha
                    , res_linewidth=res_linewidth
                    , res_linestyle=res_linestyle
                    , res_show_center=res_show_center
                    , res_markersize=res_markersize
  )

  return plt

end
# -- overloads: single record -> data
plot_eigenph!(plt, rec :: EigenphRecord; kwargs...) = plot_eigenph!(plt, rec.data; kwargs...)
plot_eigenph(rec :: EigenphRecord; kwargs...) = plot_eigenph(rec.data; kwargs...)

# -- multiple records
function plot_eigenph!(
    plt
  , records :: AbstractVector{<:EigenphRecord}
  ; cols=:all
  , selector = nothing
  , ytype :: Symbol= :phase
  , smoothwin :: Int = 1
  , ylabel = nothing
  , xlabel = "E"
  , unitsE :: Symbol = :eV
  , record_label=_default_record_label
  # -- resonance overlay control
  , idxres = nothing # nothing | :auto | Vector{IndexedResonance} | Function
  , res_detect :: NamedTuple = (;) # kwargs forwarded to find_resonances when resonances=:auto
  # -- overlay styling
  , res_color=:black
  , res_alpha :: Real=0.8
  , res_linewidth=0.8
  , res_linestyle=:solid
  , res_show_center :: Bool=false
  , res_markersize :: Real=3
  , kwargs...
  )

  # -- make sure we have stuff to plot
  isempty(records) && error("No records to plot !")

  # -- resolve columns for one record (they should be the same across records)
  data0 = records[1].data
  colidx0= resolve_cols(data0; cols=cols, selector=selector)
  isempty(colidx0) && error("No columns selected !")

  # -- decide on the resonance list
  _idxres =
    idxres === nothing ? nothing :
    idxres === :auto ? find_resonances_across_records(records; cols=colidx0, selector=nothing, res_detect...) :
    idxres isa Function ? idxres(records; cols=colidx0, selector=nothing, res_detect...) :
    idxres isa AbstractVector{<:IndexedResonance} ? idxres :
    error("idx res must be nothing | :auto | function | Vector{<:IndexedResonance}")

  # -- if we have indexed resonances, keep track of them for each record
  res_by_rec = _idxres === nothing ? nothing : resonances_by_record(_idxres, length(records))

  ylbl = ylabel === nothing ? _default_ylabel(ytype) : ylabel

  # -- loop over the records
  for (recid, r) in pairs(records)

    pfx = "$(record_label(r)) | "

    plot_eigenph!( plt, r.data
                 ; cols=colidx0
                 , selector=nothing
                 , ytype=ytype
                 , ylabel = ylbl
                 , xlabel = xlabel
                 , smoothwin=smoothwin
                 , unitsE = unitsE
                 , label=(data,col) -> pfx * data.names[col]
                 , kwargs...)

    _idxres === nothing && continue

    # -- overlay resonances from this record on the plot
    overlay_resonances!( plt, r.data, res_by_rec[recid]
                      ; ytype = ytype
                      , smoothwin=smoothwin
                      , cols_selected=colidx0
                      , unitsE = unitsE
                      , res_color=res_color
                      , res_alpha=res_alpha
                      , res_linewidth=res_linewidth
                      , res_linestyle=res_linestyle
                      , res_show_center=res_show_center
                      , res_markersize=res_markersize
                      )
  end

  return plt

end
function plot_eigenph(records :: AbstractVector{<:EigenphRecord}; kwargs...)
  plt = Plots.plot()
  return plot_eigenph!(plt, records;  kwargs...)
end

# -- overlaying resonances on eigenph data
function overlay_resonances!(
    plt
  , data :: EigenphData
  , resonances :: AbstractVector{<:Resonance}
  ; ytype ::Symbol = :phase
  , smoothwin :: Int = 1
  , unitsE :: Symbol = :eV
  , cols_selected = nothing
  , res_color = :black
  , res_alpha :: Real = 0.8
  , res_linewidth = 0.8
  , res_linestyle = :solid
  , res_show_center :: Bool = false
  , res_markersize :: Real = 3
  )

  xplot = unitsE === :au ? data.E :
          unitsE === :eV ? _convertE.(data.E, :au, :eV) :
          error("`unitsE` must be either :au or :eV")

  for r in resonances
    # -- skip non-selected columns
    (cols_selected !== nothing && !(r.col in cols_selected)) && continue

    # -- get peak position
    ys = _yseries(data, r.col; ytype=ytype, smoothwin=smoothwin)
    (ytype === :deriv || ytype === :absderiv) && (ys = _convertInvE.(ys, :au, unitsE))
    E0 = unitsE === :au ? r.E : _convertE(r.E, :au, :eV)
    y0 = linterp(xplot, ys, E0)

    # -- plot line at half max for derivatives
    (ytype === :deriv || ytype === :absderiv) && (y0/=2)

    halfwidth = ismissing(r.Γ) ? (unitsE === :eV ? 0.1 : 0.1 / AU2EV) :
                                 (unitsE === :au ? r.Γ/2 : _convertE(r.Γ/2, :au, :eV))

    # -- horizonal line indicating resonance width
    Plots.plot!( plt, [E0 - halfwidth, E0 + halfwidth], [y0, y0]
               ; label=false
               , color=res_color
               , alpha=res_alpha
               , linewidth=res_linewidth
               , linestyle=res_linestyle
              )

    res_show_center || continue

    # -- single marker indicating resonance position
    Plots.scatter!( plt, [E0], [y0]
                ; label = false
                , color=res_color
                , alpha=res_alpha
                , markersize=res_markersize
                )

  end

  return plt

end

overlay_resonances!(plt, record :: EigenphRecord, resonances; kwargs...) =
  overlay_resonances!(plt, record.data, resonances; kwargs...)

function plot_resonance!(
      plt
    , records
    , tracked
    ; key :: Symbol = :Q
    , order :: Symbol =:asc
    , connect :: Bool = false
    , include_width :: Bool =true
    , skip_missing_key :: Bool =true
    , unitsE :: Symbol = :eV # :eV or :au (Hartree)

    # -- plot aesthetics
    , xlabel = "Q"
    , ylabel = nothing
    , label = false
    , width_alpha = 0.5
    , fill_lw = 0

    # -- for adding resonance enrgies to a target state
    , targstate :: Union{Nothing, TargetState} = nothing
    , target_ref :: Symbol = :min # :min or :Q0 or :none
    , target_lw = 2
    , target_ls = :solid
    , target_label = true

    , kwargs...
    )

    nrecords = length(records)
    nrecords == length(tracked) ||
      error("Number of records ($(nrecords)) does not match the number of records in `tracked`: $(length(tracked))")

    keyvals = Float64[]
    Evals = Float64[]
    Γvals = Float64[]
    recids = Int[]

    for recid in eachindex(records)
      r = tracked[recid]
      ismissing(r) && continue

      v = meta(records[recid], key; default=missing)
      if ismissing(v)
        skip_missing_key && continue
        push!(keyvals, NaN)
      else
        push!(keyvals, float(v))
      end

      push!(Evals, r.E)
      push!(Γvals, ismissing(r.Γ) ? 0.0 : float(r.Γ))
      push!(recids, recid)
    end

    isempty(keyvals) && error("No tracked points to plot !")

    # -- sort by key
    order === :asc || order === :dsc || error("`order` must be :asc or :dsc")
    idx    = sortperm(1:length(keyvals), by = i -> keyvals[i], rev = order === :dsc)
    keyvals = keyvals[idx]
    Evals   = Evals[idx]
    Γvals   = Γvals[idx]
    recids  = recids[idx]

    # -- unit conversions, build plotting values
    halfwidth = _convertE.(0.5 .* Γvals, :au, unitsE)
    Eplot     = _convertE.(Evals, :au, unitsE)
    Γplot     = _convertE.(Γvals, :au, unitsE)

    # -- add a target state to this plot and shift resonances accordingly ?
    QV, VQ = nothing, nothing
    if targstate isa TargetState
      VQ = copy(targstate.V)
      target_ref === :min && VQ .-= minimum(VQ)
      VQ = _convertE.(VQ, :au, unitsE)
      QV = targstate.Q
      Eplot .+= [linterp(QV, VQ, q) for q in keyvals]
      targlabel = target_label === true ? targstate.name :
                  target_label === false ? false : string(target_label)
      Plots.plot!(plt, QV, VQ; label = targlabel, linestyle = target_ls, linewidth = target_lw)
    elseif targstate !== nothing
      error("targstate kwarg must be a TargetState or nothing")
    end

    # -- ylabel reflects units
    if ylabel === nothing
      ylabel = if unitsE === :eV
        "energy (eV)"
      else
        "energy (Hartree)"
      end
    end

    # -- connect the dots ?
    if connect
      Plots.plot!(plt, keyvals, Eplot; xlabel = xlabel, ylabel = ylabel, label = label, kwargs...)
      if include_width
        Elo = Eplot .- halfwidth
        Ehi = Eplot .+ halfwidth
        lc = plt.series_list[end].plotattributes[:seriescolor]
        Plots.plot!(
            plt, keyvals, Elo
          ; fillrange = Ehi
          , fillcolor = lc
          , linecolor = lc
          , fillalpha = width_alpha
          , linewidth = fill_lw
          , label = false
        )
      end
    else
      Plots.scatter!(plt, keyvals, Eplot
        ; yerror = include_width ? halfwidth : nothing
        , xlabel = xlabel, ylabel = ylabel, label = label, kwargs...
      )
    end

    return plt

end

function plot_resonance(records, tracked; kwargs...)
  plt = Plots.plot()
  return plot_resonance!(plt, records, tracked; kwargs...)
end

function plot_resonance(
      mode :: NormalMode
    , tracked
    ; targstate :: Union{TargetState, Symbol} = :auto
    , Q0 :: Real=0.0
    , kwargs...
  )
  plt = Plots.plot()
  return plot_resonance!(plt, mode.records, tracked; targstate=targstate, kwargs...)
end

function plot_resonance!(
      plt
    , mode :: NormalMode
    , tracked
    ; targstate :: Union{TargetState, Symbol} = :auto
    , Q0 :: Real = 0.0
    , kwargs...
  )
  targ = targstate === :auto ? ground_state(mode.root; Q0=Q0) :
         targstate isa TargetState ? targstate :
         nothing
  return plot_resonance!(plt, mode.records, tracked; target=targ, kwargs...)
end
