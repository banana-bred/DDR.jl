import DelimitedFiles
export resonance_table_wide
export resonance_table_long
export maketable, nrows, ncols, colindex
export widetable_to_longtable

################
#### HELPERS ###
################
@inline function putcell!(rows, i::Int, c::Int, val)
  rows[i,c] = val
  return c+1
end

###############
#### PUBLIC ###
###############

"""
    maketable(header, rows)

Makes header and rows into a Table
"""
function maketable(header :: Vector{String}, rows :: AbstractMatrix{Any})
  size(rows, 2) == length(header) ||
    error("Row dimension ($(size(rows, 2))) and header dimension ($(length(header))) do not match !")
  colmap = Dict{String, Int}(h => i for (i,h) in pairs(header))
  return Table(copy(header), Matrix{Any}(rows), colmap)
end
nrows(t::Table) = size(t.rows, 1)
ncols(t::Table) = size(t.rows, 2)
colindex(t::Table, name::AbstractString) =
  get(t.colmap, String(name), error("Could not find column '$name'"))
# -- easy indexing
Base.getindex(t::Table, i::Int, j::Int) = t.rows[i,j]
Base.getindex(t::Table, i::Int, name::AbstractString) = t.rows[i, colindex(t, name)]
# -- grab a whole column by its name
function Base.getindex(t::Table, ::Colon, name::AbstractString)
  j = colindex(t, name)
  return  t.rows[:, j]
end

function resonance_table_long(
    records :: AbstractVector{<:EigenphRecord}
  , idxres :: AbstractVector{<:IndexedResonance}
  ; key :: Symbol =:Q
  , include_geom :: Bool = true
  , include_name :: Bool = true
  , include_peak :: Bool = false
  , do_sort :: Bool = true
  , order :: Symbol = :asc
  )

  isempty(records) && error("No records to table !")
  isempty(idxres)&& error("No resonance (idxres is empty) !")

  data0 = records[1].data

  header = String["recid"]
  include_geom && push!(header, "geom")
  push!(header, string(key))
  push!(header, "col")
  include_name && push!(header, "name")
  push!(header, "E", "Γ")
  include_peak && push!(header, "peak")

  nrows = length(idxres)
  ncols = length(header)

  rows = Array{Any}(undef, nrows, ncols)

  for i in 1:nrows
    recid, r = idxres[i]
    rec = records[recid]

    c = 1
    c = putcell!(rows, i, c, recid)

    qv = meta(rec, key)
    include_geom && (c = putcell!(rows, i, c, meta(rec, :geom)))

    c = putcell!(rows, i, c, qv)
    c = putcell!(rows, i, c, r.col)

    include_name && (c = putcell!(rows, i, c, data0.names[r.col]))

    c = putcell!(rows, i, c, r.E)
    c = putcell!(rows, i, c, r.Γ)

    include_peak && (c= putcell!(rows, i, c, r.peak))
  end

    # -- optional sort
    if do_sort
        order === :asc || order === :dsc || error("order must be :asc or :desc")

        col_key = findfirst(==(string(key)), header)
        col_E   = findfirst(==("E"), header)
        col_col = findfirst(==("col"), header)

        sgn = (order === :dsc) ? -1.0 : 1.0  # keep missing last while sorting descending

        idx = sortperm(1:nrows, by = i -> begin
            q = rows[i, col_key]
            e = rows[i, col_E]
            cc = rows[i, col_col]
            recid = rows[i, 1]::Int

            qkey = ismissing(q) ? Inf : sgn * Float64(q)
            ekey = ismissing(e) ? Inf : sgn * Float64(e)

            (ismissing(q), qkey, ekey, cc, recid)
        end)

        rows = rows[idx, :]
    end

  return maketable(header, rows)

end

"""
    _homogenize_key(x; digits=4)

Homogenize floating keys (e.g., Q=0.5) for joining them across modes by rounding.
"""
@inline _homogenize_key(x :: Real; digits :: Int=4) = round(Float64(x), digits=digits)

"""
    resonance_table_wide(specs; key=:Q, order=:asc, qdigits=4, include_peak=falsr)

Build a wide table:
- row = one record (geometry, Q..)
- columns = key (+ optional geom) + for each spec: E Γ

`specs` is a vector of NamedTuples, e.g.,
  (label="BEND", records = recs, selector="singlet A1", center=0.5, window=0.05, pick=:energy, follow=true, tracked=tracked)
  `cols` can be used instead of `selector`

Returns a Table
"""
function resonance_table_wide(
    specs
  ; key :: Symbol = :Q
  , order :: Symbol = :asc
  , qdigits :: Int = 4
  , include_peak :: Bool = false
  )

  isempty(specs) && error("No specs provided for widetable !")
  order === :asc || order === :dsc || error("`order` must be :asc or :dsc")

  # -- homogenize
  resolved_specs = Vector{NamedTuple}(undef, length(specs))
  allkeys = Float64[]

  for (k, spec) in pairs(specs)
    label = String(get(spec, :label, "res$k"))
    records = get(spec, :records, nothing)
    tracked = get(spec, :tracked, nothing)

    records === nothing && error("Spec `$label` must provide `records=`")
    tracked === nothing && error("Spec `$label` must provide `tracked=`")

    length(records) == length(tracked) ||
      error("Spec '$label': length(records)=$(length(records)) must be the same as
        length(tracked)=$(length(tracked))")

    # -- collect all available keys for global row id
    for rec in records
      v = meta(rec, key; default = missing)
      ismissing(v) && continue
      push!(allkeys, _homogenize_key(v; digits=qdigits))
    end

    resolved_specs[k] = (; label, records, tracked)
  end

  isempty(allkeys) && error("Non non-missing values found for key=$key across specs !")

  # -- global grid of joined keys
  joined_keys = sort(unique(allkeys))
  order === :dsc && reverse!(joined_keys)

  # -- header
  header = String[string(key)]
  for spec in resolved_specs
    push!(header, "$(spec.label)_E")
    push!(header, "$(spec.label)_Γ")
    include_peak && push!(header, "$(spec.label)_peak")
  end

  # -- build per-spec maps: (rounded) key -> tracked resonacne OR missing
  maps = Vector{Dict{Float64, Union{Missing, Resonance}}}(undef, length(resolved_specs))
  for (k, spec) in pairs(resolved_specs)
    d = Dict{Float64, Union{Missing, Resonance}}()

    for (rec, tr) in zip(spec.records, spec.tracked)
      v = meta(rec, key; default=missing)
      ismissing(v) && continue

      qk = _homogenize_key(v; digits=qdigits)

      # -- make sure that we don't collapse duplicate rounded keys within a spec
      haskey(d, qk) && error(
        "Spec `$(spec.label)` has duplicate rounded $(key)=$(qk). " *
        "Incrase `qdigits` or take a closer look at the grids."
      )

      d[qk] = tr
    end

    maps[k] = d

  end

  # -- rows
  nrows = length(joined_keys)
  ncols = length(header)
  rows = Array{Any}(undef, nrows, ncols)

  for (irow, qk) in pairs(joined_keys)
    c = 1
    c = putcell!(rows, irow, c, qk)
    for (k, spec) in pairs(resolved_specs)
      tr = get(maps[k], qk, missing)

      if ismissing(tr)
        c = putcell!(rows, irow, c, missing)
        c = putcell!(rows, irow, c, missing)
        include_peak && (c = putcell!(rows, irow, c, missing))
        continue
      end

      c = putcell!(rows, irow, c, tr.E)
      c = putcell!(rows, irow, c, tr.Γ)
      include_peak && (c = putcell!(rows, irow, c, tr.peak))

    end
  end

  return maketable(header, rows)

end

"""
    widetable_to_longtabled(widetable, spec; drop_missing=true)

Converts a wide resonance table to a long resonance table that is more adapted to a dissociation
calcualation
"""
function widetable_to_longtable(
    widetable :: Table
  , spec :: DissCalcSpec
  ; drop_missing :: Bool = true
  )

  validate_calc_spec(widetable, spec) # from dissoc

  Qidx = _require_col(widetable, spec.qcol)

  header = String[spec.qcol, "channel", "mode", "label", "Eres", "Γ"]
  rowsv = Vector{NTuple{6, Any}}()

  for ch in spec.channels
    for cm in ch.modes
      Eidx = _require_col(widetable, "$(cm.label)_E")
      Γidx = _require_col(widetable, "$(cm.label)_Γ")

      for i in axes(widetable.rows, 1)
        Q = widetable.rows[i, Qidx]
        E = widetable.rows[i, Eidx]
        Γ = widetable.rows[i, Γidx]

        drop_missing && (ismissing(Q) || ismissing(E) || ismissing(Γ)) && continue

        push!(rowsv, (Q, ch.name, cm.mode_name, cm.label, E, Γ))
      end
    end
  end

  # -- sort by (channel, mode, Q, label)
  idx = sortperm(1:length(rowsv), by = k -> begin
    Q, ch, mode, lbl, E, Γ = rowsv[k]
    Qkey = ismissing(Q) ? Inf : Float64(Q)
    (String(ch), String(mode), Qkey, String(lbl))
  end)

  rows = Matrix{Any}(undef, length(rowsv), length(header))
  for (i_out, k_in) in enumerate(idx)
    tpl = rowsv[k_in]
    for j in 1:length(header)
      rows[i_out, j] = tpl[j]
    end
  end

  return maketable(header, rows)

end
