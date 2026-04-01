export unwrap_eigenphases!, unwrap_eigenphases
export unwrap!, unwrap
export center_records_at!
export sum_eigenphases

"""
    unwrap_eigenphases!(
        E :: AbstractVector
      , δ :: AbstractMatrix
      ; Emin = nothing
      , Emax = nothing
      , period = π
      , kwargs...
    )

Unwraps each column of δ(E) so that, at consecutive steps, the values change by at most ±period/2
within [Emin, Emax] if supplied.
Mutates δ and returns it unrwapped
"""
function unwrap_eigenphases!(
      E :: AbstractVector
    , δ :: AbstractMatrix
    ; Emin = nothing
    , Emax = nothing
    , period = π
    , kwargs...
  )

  ne = length(E)
  @assert size(δ, 1) == ne
  nδ = size(δ, 2)

  iestart = Emin === nothing ? 1 : max(1, searchsortedfirst(E, Emin))
  ieend   = Emax === nothing ? ne : min(ne, searchsortedlast(E, Emax))
  ieend - iestart <1 && return δ

  half_period = period/2

  for iδ in 1:nδ
    prev = δ[iestart, iδ]
    for ie in iestart+1 : ieend
      δ2 = δ[ie, iδ] - prev
      # -- wrap  δ2 into (-half_period, +half_period]
      δ2 = mod(δ2  + half_period, period) - half_period
      δ[ie, iδ] = prev + δ2
      prev = δ[ie, iδ]
    end
  end

  return δ

end

"""
    unwrap_eigenphases(E, δ; kwargs...) -> δ2

Unwraps each column of δ(E) so that, at consecutive steps, the values change by at most ±period/2
within [Emin, Emax] if supplied.
Does not mutate δ.
"""
unwrap_eigenphases(E, δ; kwargs...) = unwrap_eigenphases!(E, copy(δ); kwargs...)


"""
    unwrap!(data :: EigenphData; kwargs...)

Unwraps Eigenphase data.δ in place. Returns data
"""
function unwrap!(data :: EigenphData; kwargs...)
  unwrap_eigenphases!(data.E, data.δ; kwargs...)
  return data
end

"""
    unwrap!(record :: EigenphRecord; kwargs...)

Unwraps Eigenphase record.data.δ in place. Returns record
"""
function unwrap!(record :: EigenphRecord; kwargs...)
  unwrap!(record.data; kwargs...)
  return record
end

"""
    unwrap!(records :: AbstractVector{<:EigenphRecord}; kwargs...)

Unwraps Eigenphase records.data.δ in place. Returns records
"""
function unwrap!(records :: AbstractVector{<:EigenphRecord}; kwargs...)
  for r in records
    unwrap!(r; kwargs...)
  end
  return records
end

"""
    unwrap(data :: EigenphData; kwargs...)

Unwraps Eigenphase data.δ, without mutating data. Returns data
"""
function unwrap(data :: EigenphData; kwargs...)
  return EigenphData(copy(data.E), unwrap_eigenphases(data.E, data.δ; kwargs...), copy(data.names))
end

"""
    unwrap(record :: EigenphRecord; kwargs...)

Unwraps Eigenphase record.data.δ, wihtout mutating record. Returns unwrapped record
"""
function unwrap(record :: EigenphRecord; kwargs...)
  return EigenphRecord(unwrap(record.data; kwargs...), record.path, record.meta)
end

"""
    unwrap(records :: AbstractVector{<:EigenphRecord}; kwargs...)

Unwraps Eigenphase records.data.δ, wihtout mutating records. Returns unwrapped records
"""
function unwrap(records :: AbstractVector{<:EigenphRecord}; kwargs...)
  return [unwrap(r; kwargs...) for r in records]
end

"""
    linterp(E :: AbstractVector, y :: AbstractVector, x :: Real)

Linear interpolation of y(E) at x. If x is out of bounds, use the nearest endpoint.
"""
function linterp(E :: AbstractVector, y :: AbstractVector, x :: Real)

  ne = length(E)
  @assert  length(y) == ne

  # -- if x is at our outside of the bounds
  if x <= E[1]
    return y[1]
  elseif x>= E[end]
    return y[end]
  end

  ie = searchsortedlast(E, x)
  ie = clamp(ie, 1, ne-1)
  x1, x2 = E[ie], E[ie+1]
  y1, y2 = y[ie], y[ie+1]
  t = (x- x1) / (x2 - x1)

  return y1 + t*(y2 - y1)

end

"""


Shift each selected eigenphase column by adding integer multiples of the period
so that values at Eref are more or less aligned across records.

Arguments
- `cols`: `:all` or a vector of column indices
- `targtype`:
  - `:minabs` -> choose the record s.t. |δ(Eref) is smallest
  - `:median` -> chppse the record at median(δ(Eref))
"""
function center_records_at!(
    records
  ; Eref     :: Real
  , period   :: Real =π
  , cols     :: Union{Symbol,AbstractVector}=:all
  , targtype :: Symbol=:minabs
  )

  isempty(records) &&  return records

  nδ = size(records[1].data.δ, 2)
  colidx = cols === :all ? collect(1:nδ) : collect(cols)

  for j in colidx

    # -- find target δ
    vals = [linterp(r.data.E, @view(r.data.δ[:, j]), Eref) for r in records]
    if targtype == :minabs
      targ = vals[argmin(abs.(vals))]
    elseif targtype == :median
      s = sort(vals)
      targ = s[clamp(Int(ceil(length(s)/2)), 1, length(s))]
    else
      error("Unknown targtype ($(targtype))")
    end

    # -- shift by n*period
    for (r, v) in zip(records, vals)
      n = round(Int, (targ - v) / period)
      @views r.data.δ[:, j] .+= n*period
    end

  end

  return records

end

"""
    sum_eigenphases(data :: EigenphData; cols :: Union{Symbol, AbstractVector} = :all)

Return the sum of selected eigenphase columns
"""
function sum_eigenphases(data :: EigenphData; cols :: Union{Symbol, AbstractVector} = :all)
  return cols ===  :all  ? vec(sum(data.δ; dims=2)) :
                           vec(sum(@view(data.δ[:, cols]); dims=2))
end
