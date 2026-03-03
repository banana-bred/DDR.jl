"""
    fix_eigenphase_jumps(E, δ; Emin=nothing, Emax=nothing, Eref=nothing, period=π)

Fixes eigenphase (δ) jumps by ±period. If Eref is supplied, try to make all eigenphases
nearest that energy within the same branch.
"""
function _fix_eigenphase_jumps(E :: AbstractVector, δ :: AbstractArray; Emin=nothing, Emax=nothing, period=π)

  ne = length(E)
  @assert size(δ, 1) == ne
  nδ = size(δ, 2)

  iestart = Emin === nothing ? 1 : max(1,  searchsortedfirst(E, Emin))
  ieend   = Emax === nothing ? ne : min(ne, searchsortedlast(E,  Emax))

  ieend - iestart < 1 && return δ

  ahalf = period/2

  for i in 1:nδ
    prev = δ[iestart, i]
    for ie in iestart+1:ieend
      Δ = δ[ie, i] - prev
      Δ = mod(Δ+ ahalf, period) - ahalf
      δ[ie, i] = prev + Δ
      prev = δ[ie, i]
    end
  end

  return δ
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
