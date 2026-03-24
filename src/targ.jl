export pick_state, ground_state, nth_state

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
