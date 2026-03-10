import Plots
Plots.gr()

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

function plot_eigenph!(plt, data :: EigenphData
  ; cols=:all
  , selector=nothing
  , ytype :: Symbol = :phase
  , smoothwin :: Int = 1
  , label = :name
  , kwargs...
  )
  colidx = resolve_cols(data; cols=cols, selector=selector)
  isempty(colidx) && error("No columns selected")
  for c in colidx
    ys = _yseries(data, c; ytype=ytype, smoothwin=smoothwin)
    lbl = label isa Function ? label(data, c) :
      label === :name ? data.names[c] :
        label === false || label === :none || label === nothing ? false :
        string(label)
    Plots.plot!(plt, data.E,ys; label=lbl, kwargs...)
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

  # -- general plot stuff
  , ylabel = nothing
  , xlabel = "E (eV)"

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
  plt = Plots.plot(; ylabel=ylbl)

  plot_eigenph!( plt
               , data
               ; cols=cols
               , selector=nothing
               , ytype=ytype
               , smoothwin=smoothwin
               , label=label
               , kwargs...
  )

  reslist = resolve_resonances(data, resonances, colidx, res_detect)
  reslist === nothing && return plt

  overlay_resonances!(plt, data, reslist
                    ; ytype=ytype
                    , smoothwin=smoothwin
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
function plot_eigenph(
    records :: AbstractVector{<:EigenphRecord}
  ; cols=:all
  , selector = nothing
  , ytype :: Symbol= :phase
  , smoothwin :: Int = 1
  , ylabel = nothing
  , xlabel = "E (eV)"
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
  plt = Plots.plot(; ylabel=ylbl)

  # -- loop over the records
  for (recid, r) in pairs(records)

    pfx = "$(record_label(r)) | "

    plot_eigenph!( plt, r.data
                 ; cols=colidx0
                 , selector=nothing
                 , ytype=ytype
                 , smoothwin=smoothwin
                 , label=(data,col) -> pfx * data.names[col]
                 , kwargs...)

    _idxres === nothing && continue

    # -- overlay resonances from this record on the plot
    overlay_resonances!( plt, r.data, res_by_rec[recid]
                      ; ytype = ytype
                      , smoothwin=smoothwin
                      , cols_selected=colidx0
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

# -- overlaying resonances on eigenph data
function overlay_resonances!(
    plt
  , data :: EigenphData
  , resonances :: AbstractVector{<:Resonance}
  ; ytype ::Symbol = :phase
  , smoothwin :: Int = 1
  , cols_selected = nothing
  , res_color = :black
  , res_alpha :: Real = 0.8
  , res_linewidth = 0.8
  , res_linestyle = :solid
  , res_show_center :: Bool = false
  , res_markersize :: Real = 3
  )

  for r in resonances
    # -- skip non-selected columns
    (cols_selected !== nothing && !(r.col in cols_selected)) && continue

    # -- get peak position
    ys = _yseries(data, r.col; ytype=ytype, smoothwin=smoothwin)
    y0 = linterp(data.E, ys, r.E)

    # -- plot line at half max for derivatives
    (ytype === :deriv || ytype === :absderiv) && (y0/=2)

    halfwidth = ismissing(r.Γ) ? 0.1 : r.Γ/2
    # halfwidth = r.Γ/2

    # -- horizonal line indicating resonance width
    Plots.plot!( plt, [r.E - halfwidth, r.E + halfwidth], [y0, y0]
               ; label=false
               , color=res_color
               , alpha=res_alpha
               , linewidth=res_linewidth
               , linestyle=res_linestyle
              )

    res_show_center || continue

    # -- single marker indicating resonance position
    Plots.scatter!( plt, [r.E], [y0]
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

    # -- plot aesthetics
    , xlabel = "Q"
    , ylabel = "E (eV)"
    , label = false
    , width_alpha = 0.5
    , fill_lw = 0

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
    perm    = sortperm(1:length(keyvals), by = i -> keyvals[i], rev = order === :dsc)
    keyvals = keyvals[perm]
    Evals   = Evals[perm]
    Γvals   = Γvals[perm]
    recids  = recids[perm]

    halfwidth = 0.5 .* Γvals

    if connect
      Plots.plot!(plt, keyvals, Evals; xlabel = xlabel, ylabel = ylabel, label = label, kwargs...)
      if include_width
        Elo = Evals .- halfwidth
        Ehi = Evals .+ halfwidth
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
      Plots.scatter!(plt, keyvals, Evals
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
