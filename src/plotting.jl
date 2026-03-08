import Plots
Plots.gr()

export plot_eigenph, plot_eigenph!
export overlay_resonances!

#################
#### HELPERS ####
#################

"""
    _yseries(data,  col; ytype, smoothwin)

returns the appropriate ytype series given Eigenphase data
- ytype :: Symbol
  - :phase (δ)
  - :deriv (dδ/dE)
  - :absderiv (|dδ/dE|)
- smoothwin :: Int -> moving average smoothing window
"""
function _yseries(
    data      :: EigenphData
  , col       :: Int
  ; ytype     :: Symbol = :phase
  , smoothwin :: Int  = 1
  )
  δ = @view data.δ[:,col]
  ytype === :phase && return δ
  dδ = _dydE(data.E, δ)
  dδ = smoothwin > 1 ? _moving_average(dδ, smoothwin) : dδ
  ytype === :deriv && return dδ
  ytype === :absderiv && return abs.(dδ)
  error("Symbol ytype must be :phase, :deriv, or :absderiv")
end

function _default_record_label(record :: EigenphRecord)
  g = geom_index(record.path)
  return g === nothing  ? splitpath(record.path)[end] : "geom = $(g)"
end

function _default_ylabel(ytype :: Symbol)
  ytype === :phase    && return "δ"
  ytype === :deriv    && return "dδ/dE"
  ytype === :absderiv && return "|dδ/dE|"
  error("Symbol ytype must be :phase, :deriv, or :absderiv")
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
    x
  ; resonances=nothing
  , cols=:all
  , selector=nothing
  , ytype :: Symbol = :phase
  , smoothwin :: Int = 1
  , ylabel=nothing
  , res_tick_width=nothing
  , res_color=:black
  , res_alpha :: Real=0.8
  , res_linewidth=0.8
  , res_linestyle=:solid
  , res_show_center :: Bool=true
  , res_markersize :: Real=3
  , kwargs...
  )
  data = x isa EigenphRecord ? x.data : x
  ylbl = ylabel === nothing ? _default_ylabel(ytype) : ylabel
  plt = Plots.plot(; ylabel=ylbl)
  plot_eigenph!(plt, data;ytype=ytype, smoothwin=smoothwin, cols=cols, selector=selector, kwargs...)
  resonances === nothing && return plt
  colidx = resolve_cols(data; cols=cols, selector=selector)
  overlay_resonances!(plt, data, resonances
                    ; ytype=ytype
                    , smoothwin=smoothwin
                    , cols_selected=colidx
                    , res_tick_width=res_tick_width
                    , res_color=res_color
                    , res_alpha=res_alpha
                    , res_linewidth=res_linewidth
                    , res_linestyle=res_linestyle
                    , res_show_center=res_show_center
                    , res_markersize=res_markersize
                    )
  return plt
end

# -- overloads for a single record
plot_eigenph!(plt, rec :: EigenphRecord; kwargs...) = plot_eigenph!(plt, rec.data; kwargs...)
plot_eigenph(rec :: EigenphRecord; kwargs...) = plot_eigenph(rec.data; kwargs...)

# -- multiple records
function plot_eigenph(
    records :: AbstractVector{<:EigenphRecord}
  ; record_label=_default_record_label
  , idxres = nothing #  nothing | :auto | Vector{IndexedResonance} | Function
  , res_detect :: NamedTuple = (;) # forwarded to find_resonances_across_records
  # -- what to plot
  , ytype:: Symbol= :phase
  , smoothwin :: Int = 1
  , ylabel = nothing
  , cols = :all
  , selector = nothing
  # -- resonance overlay styling
  , res_tick_width = nothing
  , res_color = :black
  , res_alpha :: Real = 0.8
  , res_linewidth = 0.8
  , res_linestyle = :solid
  , res_show_center :: Bool = true
  , res_markersize :: Real = 3
  , kwargs...
  )

  # -- make sure we have stuff to plot
  isempty(records) && error("No records to plot !")

  # -- resolve columns for one record (they should be the same across records)
  data0 = records[1].data
  colidx0= resolve_cols(data0; cols=cols, selector=selector)
  isempty(colidx0) && error("No columns selected !")

  # -- decide on the resonance list
  _idxres = idxres
  if _idxres === :auto
    # -- automatically detect resonances per record, using the same columns as the plot,
    _idxres = find_resonances_across_records(records;cols=colidx0, selector=nothing, res_detect...)
  elseif _idxres isa Function
    _idxres = _idxres(records; cols=colidx0, selector=nothing, res_detect...)
  elseif !(_idxres === nothing || _idxres isa AbstractVector{<:IndexedResonance})
    error("idxres must be nothing, :auto, a Function, or a vector of IndexedResonance")
  end

  # -- if we have indexed resonances, keep track of them for each record
  res_by_rec = nothing
  if _idxres !== nothing
    res_by_rec = [Resonance[] for _ in 1:length(records)]
    for (recid, res) in _idxres
      push!(res_by_rec[recid], res)
    end
  end

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
                      , res_tick_width=res_tick_width
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
  , res_tick_width = nothing
  , res_color = :black
  , res_alpha :: Real = 0.8
  , res_linewidth = 0.8
  , res_linestyle = :solid
  , res_show_center :: Bool = true
  , res_markersize :: Real = 3
  )

  tw = res_tick_width === nothing ? 4*median(diff(data.E)) : res_tick_width

  for r in resonances
    # -- skip non-selected columns
    (cols_selected !== nothing && !(r.col in cols_selected)) && continue

    # -- get peak position
    ys = _yseries(data, r.col; ytype=ytype, smoothwin=smoothwin)
    y0 = linterp(data.E, ys, r.E)

    # -- plot line at half max for derivatives
    ytype === :deriv || ytype === :absderiv && y0/=2

    halfwidth = ismissing(r.Γ) ? tw/2 : r.Γ/2
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
