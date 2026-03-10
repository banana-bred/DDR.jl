import DelimitedFiles
export show_table
export write_table
export resonance_table_wide
export resonance_table_long

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

# -- show a table
function show_table(header, rows; maxrows=30)
    n = min(size(rows,1), maxrows)
    println(join(header, "  |  "))
    println("-"^min(120, 6*length(header)+40))
    for i in 1:n
        println(join(string.(rows[i, :]), "  |  "))
    end
    if size(rows,1) > n
        println("... $(size(rows,1)-n) more")
    end
end

"""
    resonance_table_wide(records idxres; key=:Q, sortby=key, order=:asc
      , specs, include_geom=true)

Build a wide table:
- row = one record (geometry, Q..)
- columns = key (+ optional geom) + for each spec: E Γ

`specs` is a vector of NamedTuples, e.g.,
  (label="BEND", selector="singlet A1", center=0.5, window=0.05, pick=:maxpeak, follow=true)
  `cols` can be used instead of `selector`

Returns (header :: Vector{String}, rows :: Matrix{Any})
"""
function resonance_table_wide(
    records :: AbstractVector{<:EigenphRecord}
  , idxres :: AbstractVector{<:IndexedResonance}
  ; key :: Symbol = :Q
  , sortby :: Symbol = key
  , order :: Symbol = :asc
  , specs
  , include_geom :: Bool = true
  , include_peak :: Bool = false
  )

  isempty(records) && error("No records to table !")
  nrecords = length(records)
  data0 = records[1].data

  # -- group resonances together
  byrec = resonances_by_record(idxres, length(records))

  # -- record ordering
  keyvals = [meta(r, sortby) for r in records]
  ord = sortperm(1:nrecords, by = i ->( ismissing(keyvals[i]), keyvals[i] ), rev = (order === :desc))

  # -- resolve spec -> (label, col, tracked vector)
  resolved = Vector{NamedTuple}(undef, length(specs))
  for (k, spec) in pairs(specs)
    label = get(spec, :label, "res $k")

    col =
      haskey(spec, :col) ? spec.col :
        haskey(spec, :selector) ? begin
          cols = resolve_cols(data0; cols=:all, selector=spec.selector)
          length(cols) == 1 || error("Spec '$label' selector must match exactly 1 column; got $(length(cols))")
          cols[1]
        end :
        error("Spec '$label' needs :col or :selector")

    center = haskey(spec, :center) ? spec.center : error("$Spec '$label' needs :center")
    window = haskey(spec, :window) ? spec.window : error("$Spec '$label' needs :window")
    follow = haskey(spec, :follow) ? spec.follow : true
    follow_window = haskey(spec, :follow_window) ? spec.follow_window : window
    pick = haskey(spec, :pick) ? spec.pick : :closest

    tr = track_resonance( records
                        , byrec
                        ; col=col
                        , center=center
                        , window=window
                        , follow=follow
                        , follow_window=follow_window
                        , pick=pick
    )
    resolved[k] = (; label,col, tr)
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

  # -- rows
  nrows = length(records)
  ncols = length(header)
  rows = Array{Any}(undef, nrows, ncols)

  for (irow, recid) in pairs(ord)
    rec = records[recid]
    c = 1

    rows[irow,  c] = meta(rec, key) ; c += 1

    if include_geom
      rows[irow, c] = meta(rec, :geom) ; c += 1
    end

    for spec in resolved
      rsel = spec.tr[recid]
      if ismissing(rsel)
        rows[irow, c] = missing ; c += 1
        rows[irow, c] = missing ; c += 1
        if include_peak
          rows[irow, c] = missing ; c += 1
        end
      else
        rows[irow, c] = rsel.E ; c += 1
        rows[irow, c] = rsel.Γ ; c += 1
        if include_peak
          rows[irow, c] = rsel.peak ; c += 1
        end
      end
    end

  end

  return header, rows

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

  return header, rows

end

"""
    write_table(path, header, rows; delim=',')

Writes a header row + `rows` (matrix | vector of vectors) to disk.
Returns `path`.
"""
function write_table(
    path  :: AbstractString
  , header :: Vector{String}
  , rows
  ; delim=','
  )
  open(path, "w") do io
    DelimitedFiles.writedlm(io, permutedims(header), delim)
    DelimitedFiles.writedlm(io, rows, delim)
  end
  return path
end
