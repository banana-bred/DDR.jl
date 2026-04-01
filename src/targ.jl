using QuadGK: quadgk

export pick_state, ground_state, nth_state
export ζ_gauss, ζ_path_fit, ζ_Q
export smin_V, E0_harmonic

function pick_state(targE :: TargetEnergies, col :: Int)
  1 <= col <= size(targE.E, 2) || error("col=$col out of range [1:$(size(targE.E, 2))]")
  return TargetState(copy(targE.Q), vec(@view targE.E[:, col]), targE.names[col], col)
end

"""
    ground_state(targE; Q0=0.0)

Pick the lowest energy state at the grid point closes to Q0,
then return that column as a TargetState struct
"""
function ground_state(targE :: TargetEnergies; Q0 :: Real = 0.0)
  i0 = argmin(abs.(targE.Q .- float(Q0)))
  j0 = argmin(@view(targE.E[i0, :]))
  return pick_state(targE, j0)
end

"""
    ground_state(mode_root, Q0 = 0.0, fname = "target.energies")

Wrapper for `ground_state`. Reads the file then returns a TargetState for the ground state
"""
function ground_state(mode_root :: AbstractString; Q0 :: Real = 0.0, fname :: AbstractString = "target.energies")
  path = joinpath(mode_root, "target.energies")
  isfile(path) || return nothing
  targE = read_target_energies(path)
  return ground_state(targE; Q0=Q0)
end

"""
    nth_state(targE; n=1, Q0=0.0)

Pick the state with the given n at the grid point closest to Q0.
- n=1 -> ground (ground by this definition)
- n=2 -> first excited (at Q0)
Returns that column as a TargetState.
"""
function nth_state(targE :: TargetEnergies; n :: Int = 1, Q0 :: Real = 0.0)
    n >= 1 || error("n must be ≥ 1")
    i0 = argmin(abs.(targE.Q .- float(Q0)))
    e0 = collect(@view targE.E[i0, :])
    n <= length(e0) || error("n=$n > nstates=$(length(e0))")

    # sort column indices by energy at Q0
    ord = sortperm(e0)
    j = ord[n]
    return pick_state(targE, j)
end

"""
    ζ_gauss(pathfit, s; α=1.0, s0=nothing)

Returns gaussian for the ground state wavefunction.
If s0 is supplied, use that as the minimum.
Otherwise, try and determine it from V(s)

    exp(-α * (s-s0)² / 2)
"""
function ζ_gauss(pathfit, s; α=1.0, s0 :: Union{Real, Nothing} = nothing)
  s0 = s0 === nothing ? smin_V(pathfit) : Float64(s0)
  x = Float64(s) - s0
  return exp(-0.5 * Float64(α) * x^2)
end

"""
    ζ_Q(target, Q; Q0=nothing)

Normalized harmonic ground-state wavefunction for one mode in dimensionless Q.
If Q0 is not supplied, it is taken as the minimum of the target potential.
"""
function ζ_Q(target :: TargetState, Q :: Real; Q0 :: Union{Real, Nothing} = nothing)
  q0 = Q0 === nothing ? target.Q[argmin(target.V)] : Float64(Q0)
  x = Float64(Q) - q0
  return exp(-0.5*x^2) / π^0.25
end

"""
    ζ_path_fit(pathfit, targets_by_mode; kind=:linear)

Builds a normalized path wavefunction ζ(s) from the product of single-mode wavefunctions
evaluated along the sampled path, then interpolated in s
"""
function ζ_path_fit(
    pathfit :: DissociationPathwayFit
  , targets_by_mode :: AbstractDict{String, <:TargetState}
  ; kind :: Symbol = :linear
  )

  modes = pathfit.path.modes
  sgrid = pathfit.s
  qmat  = pathfit.q

  ns = length(sgrid)
  vals = Vector{Float64}(undef, ns)

  for is in 1:ns
    amp = 1.0
    for (imode, mode_name) in enumerate(modes)
      haskey(targets_by_mode, mode_name) || error("No target state found for mode '$mode_name'")
      amp *= ζ_Q(targets_by_mode[mode_name], qmat[is, imode])
    end
    vals[is] = amp
  end

  raw_fit = _make_interp(sgrid, vals; kind=kind, clampx=false)

  smin, smax = extrema(sgrid)
  norm2, _ = quadgk(s -> abs2(raw_fit(s)), smin, smax)
  norm2 > 0 || error("Path wavefunction has zero norm on the sample pathway")

  return s -> raw_fit(s) / sqrt(norm2)

end

"""
    smin_V(pathfit :: DissociationPathwayFit)

Approximates the minimal V(s) on the fitted dissociation path
"""
function smin_V(pathfit :: DissociationPathwayFit)
  vals = pathfit.V_fit.(pathfit.s)
  return pathfit.s[argmin(vals)]
end

"Harmonic oscillator energy"
E0_harmonic(ω :: Real, n :: Int = 0) = (n + 0.5) * Float64(ω)
