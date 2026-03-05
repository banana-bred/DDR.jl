import Plots
Plots.gr()

export plot_eigenph, plot_records

"""
    function plot_eigenph(data; cols=:all, selector=nothing,  label=:name, ylabel="δ". kwargs...) -> Plot

Plot eigenphases δ(E) from a given EigenphData

-  Choose columns by `cols` or `selector`
- `label=:name` uses data.names[col] as the legend label
"""
function plot_eigenph(
      data
    ; cols=:all
    , selector=nothing
    , label=:name
    , ylabel :: AbstractString = "δ"
    , kwargs...
  )
  colidx = resolve_cols(data; cols=cols, selector=selector)
  isempty(colidx) && error("No columns selected")
  plt = Plots.plot(; kwargs...)
  for c in colidx
    lbl = label === :name ? data.names[c] :
        label === false || label === :none ? false : string(label)
        Plots.plot!(plt,  data.E, @view(data.δ[:, c]); label=lbl)
  end
  return plt
end

function _default_record_label(record :: EigenphRecord)
  g = geom_index(record.path)
  return g === nothing  ? splitpath(record.path)[end] : "geom = $(g)"
end

"""
    function plot_records(
          records :: AbstractVector{<:EigenphRecord}
        ; cols=:all
        , selector=nothing
        , record_label=_default_record_label
        , kwargs...
      )

Plots eigenphases for multiple records.
Useful for plotting the same data across geometries.

- `selector`/`cols` select eigenphase columns
- `record_label` is a function (record -> String) that tells us how to label data
"""
function plot_records(
      records :: AbstractVector{<:EigenphRecord}
    ; cols=:all
    , selector=nothing
    , record_label=_default_record_label
    , kwargs...
  )
  isempty(records) && error("No records to plot !")
  data = records[1].data
  colidx = resolve_cols(data; cols=cols, selector=selector)
  isempty(colidx) && error("No columns selected !")
  plt = Plots.plot(; kwargs...)
  for r in records
    for c in colidx
      lbl = "$(record_label(r)) | $(r.data.names[c])"
      Plots.plot!(plt, r.data.E, @view(r.data.δ[:, c]); label=lbl)
    end
  end
  return plt
end
