import DelimitedFiles
export resonance_table_wide
export resonance_table_long
export maketable, nrows, ncols, colindex

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

"""
    resonance_table_wide(records idxres; key=:Q, sortby=key, order=:asc
      , specs, include_geom=true)

Build a wide table:
- row = one record (geometry, Q..)
- columns = key (+ optional geom) + for each spec: E Γ

`specs` is a vector of NamedTuples, e.g.,
  (label="BEND", selector="singlet A1", center=0.5, window=0.05, pick=:energy, follow=true, tracked=tracked)
  `cols` can be used instead of `selector`

Returns a Table
"""
function resonance_table_wide(
    records :: AbstractVector{<:EigenphRecord}
  , specs
  ; key :: Symbol = :Q
  , sortby :: Symbol = key
  , order :: Symbol = :asc
  , include_geom :: Bool = true
  , include_peak :: Bool = false
  , do_sort :: Bool = true
  , kwargs...
  )

  isempty(records) && error("No records to table !")
  nrecords = length(records)

  # -- resolve spec -> (label, col, tracked vector)
  resolved = Vector{NamedTuple}(undef, length(specs))
  for (k, spec) in pairs(specs)

    label = get(spec, :label, "res $k")
    tr = get(spec, :tracked, nothing)
    tr === nothing && error("Spec '$label' must provide `tracked=` (Vector{Union{Resonance, Missing}})")

    length(tr) == nrecords || error("Spec '$label' has $(length(tr)) records, but $nrecords are expected")

    resolved[k] = (; label, tr)

  end

  # -- header
  header = String[]
  push!(header, string(key))
  include_geom && push!(header, "geom")
  for spec in resolved
    push!(header, "$(spec.label)_E")
    push!(header, "$(spec.label)_Γ")
    include_peak && push!(header, "$(spec.label)_peak")
  end

  # -- record ordering
  ord = collect(1:nrecords)
  if do_sort
    (order === :asc || order === :dsc) || error("`order` must be :asc or :dsc")
    sgn = (order === :dsc) ? -1.0 : 1.0
    ord = sortperm(1:nrecords, by = i ->
      begin
        v = meta(records[i], sortby; default=missing)
        vkey = ismissing(v) ? Inf : sgn * Float64(v)
        (ismissing(v), vkey)
      end)
  end

  # -- rows
  nrows = length(records)
  ncols = length(header)
  rows = Array{Any}(undef, nrows, ncols)

  for (irow, recid) in pairs(ord)
    rec = records[recid]
    c = 1

    c = putcell!(rows, irow, c, meta(rec, key))

    include_geom && (c = putcell!(rows, irow, c, meta(rec, :geom)))

    for spec in resolved
      rsel = spec.tr[recid]
      if ismissing(rsel)
        c = putcell!(rows, irow, c, missing)
        c = putcell!(rows, irow, c, missing)
        include_peak && (c = putcell!(rows, irow, c, missing))
      else
        c = putcell!(rows, irow, c, rsel.E)
        c = putcell!(rows, irow, c, rsel.Γ)
        include_peak && (c = putcell!(rows, irow, c, rsel.peak))
      end
    end

  end

  return maketable(header, rows)

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

        perm = sortperm(1:nrows, by = i -> begin
            q = rows[i, col_key]
            e = rows[i, col_E]
            cc = rows[i, col_col]
            recid = rows[i, 1]::Int

            qkey = ismissing(q) ? Inf : sgn * Float64(q)
            ekey = ismissing(e) ? Inf : sgn * Float64(e)

            (ismissing(q), qkey, ekey, cc, recid)
        end)

        rows = rows[perm, :]
    end

  return maketable(header, rows)

end
