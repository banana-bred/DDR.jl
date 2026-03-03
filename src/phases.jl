export unwrap_eigenphases!
export unwrap_records!
export center_records_at
export eigenphase_sum

"""
    unwrap_eigenphases!(
        E :: AbstractVector
      , δ :: AbstractVector
      ; Emin = nothing
      , Emax = nothing
      , period = π
    )

Unwraps each column of δ(E) so that, at consecutive steps, the values change by at most ±period/2
within [Emin, Emax] if supplied.
Mutates δ and returns it unrwapped
"""
function unwrap_eigephases!(
      E :: AbstractVector
    , δ :: AbstractVector
    ; Emin = nothing
    , Emax = nothing
    , period = π
  )

  ne = length(E)
  @assert size(δ, 1) == ne
  nδ = size(δ, 2)

  iestart = Emin === nothing ? 1 : max(1, searchsortedfirst(E, Emin))
  ieend   = Emax === nothing ? ne : min(ne, searchsortedfirst(E, Emax))
  ieend - iestart <1 && return δ

  half_period = period/2

  for  i in 1:nδ
    prev = δ[iestart, j]
    for ie in iestart+1 : ieend
      δ2 = δ[i, ie] - prev
      # -- wrap  δ2 into (-half_period, +half_period]
      δ2 = mod(δ2  + half_period, period) - half_period
      δ[i, ie] = prev + δ2
      prev = δ[i, ie]
    end
  end

  return δ

end

"""
    unwrap_records!(records; kwargs...)

Apply `unwrap_eigenphases!` to each record, in place.
"""
@inline function unwrap_records!(records; kwargs...)
    for r in records
        unwrap_eigenphases!(r.data.E, r.data.δ; kwargs...)
    end
    return records
end

"""
    linterp_at(E :: AbstractVector, y :: AbstractVector, x :: Real)

Linear interpolation of y(E) at x. If x is out of bounds, use the nearest endpoint.
"""
function  linterp_at(E :: AbstractVector, y :: AbstractVector, x :: Real)

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
  x1, x2 = E[i], E[i+1]
  y1, y2 = y[i], y[i+1]
  t = (x- x1) / (x2 - x1)

  return y1 + t*(y2 = y1)

end

@@@
function center_records_at!(
    records
  ; Eref
  , period=π
  , cols=:all
  , choose_targ=:minabs
  )
end

# # -- centering eigenphase
# Eref = .01
# for i=2:size(data[1], 2)
#   E = data[1][:,1]
#   ieref = argmin(abs.(E.-Eref))
#   vals = [data[g][ieref, i] for g in 1:ngeoms]
#   target = minimum(abs.(vals))
#   for igeom in 1:ngeoms
#     n = round(Int, (target - data[igeom][ieref, i]) / period)
#     data[igeom][:, i] .+= n*period
#   end
# end
