using QuadGK: quadgk

export find_s_all_Eres, find_s_all_U, find_s_all
export find_s_Eres, find_s_U, find_s
export select_root
export thermal_rate, thermal_rates
export σ
export channel_E0
export get_spinmult
export resolve_spinwght

###############
### HELPERS ###
###############

@inline _family_name(name :: AbstractString) = String(replace(name, r"\s*\[\d+\]$" => ""))

@inline _branch_from_s(s::Real) =
    s < 0 ? :negative :
    s > 0 ? :positive :
            :inner

"""
    get_spinmult(label)

Get the spin multiplicity from a label, e.g., 'singlet' -> 1, 'triplet' -> 3.
If this doesn't work, try the form, e.g., 1A1 -> 1, 3B2 -> 3.
Else, return nothing
"""
function get_spinmult(label :: AbstractString)
  s = lowercase(strip(String(label)))
  # -- words
  occursin("singlet", s) && return 1
  occursin("doublet", s) && return 2
  occursin("triplet", s) && return 3
  occursin("quartet", s) && return 4
  occursin("quintet", s) && return 5
  occursin("sextet", s) && return 6
  # -- numeric
  m = match(r"^\s*(\d+)\s*[A-Za-z]", String(label))
  return m === nothing ? nothing : parse(Int, m.captures[1])
end

"""
    resolve_spinwght(family, spinmult, target_spinmult, spin_mode)

Resolves the resonance spin multiplicity and statistical weight for an
unpolarized electron incident on a target of multiplicity `target_spinmult`

`spin_mode`
----------
- :auto   -> infer `M` from `family` unless `spinmult` is supplied
- :manual -> use `spinmult` directly
- :none   -> disable spin weighting, return `(nothing, 1.0)`

Returns `(M,w)`, where
- `M`: resonance multiplicity 2S+1
- `w = M/(2*target_spinmult)`

If `target_spinmult` is not give, returns weigh `1.0`.
"""
function resolve_spinwght(
      family          :: AbstractString
    , spinmult        :: Union{Int, Nothing} = nothing
    , target_spinmult :: Union{Int, Nothing} = nothing
    , spin_mode       :: Symbol = :auto
  )

  M = if spin_mode === :none
    nothing
  elseif spin_mode === :manual
    spinmult
  elseif spin_mode === :auto
    isnothing(spinmult) ? get_spinmult(family) : spinmult
  else
    throw(ArgumentError("unsupported spin_mode=$spin_mode"))
  end

  isnothing(M) && return (nothing, 1.0)
  isnothing(target_spinmult) && return (M, 1.0)

  # -- unpolarized electron
  w = Float64(M) / (2.0 * Float64(target_spinmult))
  return (M, w)

end

####################
### ROOT FINDERS ###
####################

"""
    find_s_all(pathfit, f, target; tol=1e-10)

Solves `f(pathfit, s) = target`.

Fins all roots all corresponding `s`
"""
function find_s_all(
    pathfit :: DissociationPathwayFit
  , f
  , target :: Float64
  ; tol :: Float64 = 1e-10
  )

  y = [f(pathfit, s) for s in pathfit.s]
  br = _find_brackets(pathfit.s, y, target)
  isempty(br) && return Float64[]
  # isempty(br) && error("No solution for requested target value ($target) on the fitted path")

  roots = Float64[]

  for (i, j) in br
    if i == j
      push!(roots, Float64(pathfit.s[i]))
    else
      g(s) = f(pathfit, s) - target
      push!(roots, _bisect_root(g, pathfit.s[i], pathfit.s[j]; tol=tol))
    end
  end

  sort!(roots)
  unique!(roots)
  return roots

end
@inline find_s_all_Eres(pathfit :: DissociationPathwayFit, ε :: Real; kwargs...) =
  find_s_all(pathfit, Eres, Float64(ε); kwargs...)
@inline find_s_all_U(pathfit :: DissociationPathwayFit, Etot :: Real; kwargs...) =
  find_s_all(pathfit, U, Float64(Etot); kwargs...)

"""
    find_s(pathfit, f, target; branch=:inner, sref=nothing, tol=1e-10)

Solve `f(pathfit, s) = target` and picks the best `s` according to branch:
- :first    -> picks the first root
- :last     -> picks the last root
- :inner    -> picks the root with smallest magnitude (i.e. inner; closest to `s=0`)
- :negative -> picks the negative root with smallest magnitude
- :positive -> picks the positive root with smallest magnitude
- :closest  -> picks the root that is closest to sref, which must be supplied in this case
"""
function find_s(
      pathfit :: DissociationPathwayFit
    , f
    , target :: Float64
    ; branch :: Symbol = :inner
    , sref :: Union{Nothing, Real} = nothing
    , tol :: Float64 = 1e-10
  )
  roots = find_s_all(pathfit, f, target; tol=tol)
  isempty(roots) && error("No solution for requested target value $target on the fitted path")
  s = select_root(roots; branch=branch, sref=sref)
  s === nothing && error("None of the found roots matched for branch=$branch for target=$target")
  return s
end
@inline find_s(pathfit::DissociationPathwayFit, f, target::Real; kwargs...) =
  find_s(pathfit, f, Float64(target); kwargs...)
@inline find_s_Eres(pathfit::DissociationPathwayFit, ε::Real; kwargs...) =
  find_s(pathfit, Eres, ε; kwargs...)
@inline find_s_U(pathfit::DissociationPathwayFit, Etot::Real; kwargs...) =
  find_s(pathfit, U, Etot; kwargs...)

"""
    select_root(roots; branch, sref)

Pick a root from `roots` based on `branch`:
- :first    -> picks the first root
- :last     -> picks the last root
- :inner    -> picks the root with smallest magnitude (i.e. inner; closest to `s=0`)
- :negative -> picks the negative root with smallest magnitude
- :positive -> picks the positive root with smallest magnitude
- :closest  -> picks the root that is closest to sref, which must be supplied in this case
"""
function select_root(
    roots :: AbstractVector{<:Real}
  ; branch :: Symbol = :inner
  , sref :: Union{Nothing, Real} = nothing
  )
  isempty(roots) && return nothing
  r = Float64.(roots)
  if branch === :first
    return first(r)
  elseif branch === :last
    return last(r)
  elseif branch === :inner
    return r[argmin(abs.(r))]
  elseif branch === :negative
    neg = filter(s -> s < 0.0, r)
    isempty(neg) && return nothing
    return maximum(neg) # closest to 0
  elseif branch === :positive
    pos = filter(s -> s > 0.0, r)
    isempty(pos) && return nothing
    return minimum(pos) # closest to 0
  elseif branch === :closest
    sref === nothing && throw(ArgumentError("branch=:closest requires non-nothing sref"))
    return r[argmin(abs.(r .- Float64(sref)))]
  else
    throw(ArgumentError("unsupported branch choice $branch"))
  end

end

#########################
### THERMAL AVERAGING ###
#########################

"""
    thermal_rate(energies, σ, T; interp_kind=:linear, m=1.0)

Thermal (state-selected) average over electron speeds at a given temperature T, assuming a Maxwellian distribution

    k(T) = sqrt(8/(π*m)) * (T)^(-3/2) * ∫ σ(E) E exp(-E/T) dE (atomic units)

- m: electron mass (1.0 in atomic units)
- Input:  atomic units
- Output: atomic units
"""
function thermal_rate(
    energies :: AbstractVector{<:Real}
  , σ :: AbstractVector{<:Real}
  , T :: Real
  ; interp_kind :: Symbol = :linear
  , m :: Real = 1.0 # electron mass
  )

  nE, nσ = length.((energies, σ))
  nE == nσ || error("energies/σ length mismatch ($nE/$nσ)")

  T > 0 || error("Temperature ($T) must be positive")

  T = Float64(T)
  E_, σ_ = Float64.(energies), Float64.(σ)

  σfit = _make_interp(E_, σ_; kind=interp_kind)
  Emin, Emax = extrema(E_)
  prefactor = sqrt(8.0/(π*m)) * T^(-1.5)

  integrand(E) = σfit(E) * E * exp(-E/T)
  vals, _ = quadgk(integrand, Emin, Emax)

  return prefactor * vals

end

"""
    thermal_average(energies, σ, Tgrid; kwargs...)

Determine state-selected thermal rate coefficients on a temperature grid.
"""
function thermal_rates(
    energies :: AbstractVector{<:Real}
  , σ :: AbstractVector{<:Real}
  , Tgrid :: AbstractVector{<:Real}
  ; kwargs...
  )
  return [thermal_rate(energies, σ, T; kwargs...) for T in Tgrid]
end

######################
### CROSS SECTIONS ###
######################

"""
    channel_E0(channel, mode_freqs)

Determine the harmonic zero-point energy E0 from a subset of mode frequencies
"""
function channel_E0(channel :: Channel, mode_freqs :: AbstractDict{String, <:Real})
  return 0.5 * sum(Float64(mode_freqs[cm.mode_name]) for cm in channel.modes)
end

"""
    σ(pathfit, ε; E0, ζ_eval, Etot=nothing, branch_ε=:inner, branch_Etot=:same_side, tol=1e-10)

Compute the dissociation cross section for a single electron energy `ε`.

Arguments
---------
- `pathfit`: the PathFit for this dissociative pathway
- `ε`: electron energy
- `E0`: the combined zero-point energy of the contributing modes.
- `ζ_eval`: the ground-state wavefunction evaluator
- `Etot`: the total energy. If nothing, this is just `E0 + ε`
- `branch_ε`: how to choose s_ε. See `find_s` for options here
- `branch_Etot`: very similar to `branch_ε`, but for finding s_E. Same options as `branch_ε` exist, in addition to
  - `:same_side`:  picks the same branch as `s_ε`: results in `:negative`, `:positive`, or `:inner`
  - `:closest_to_sε`: picks the root closest to s_ε
"""
function σ(
    pathfit :: DissociationPathwayFit
  , ε :: Real
  ; E0 :: Real
  , ζ_eval
  , Etot :: Union{Real, Nothing} = nothing
  , branch_ε :: Symbol = :inner
  , branch_Etot :: Symbol = :same_side
  , tol :: Real = 1e-10
  )

  ε = Float64(ε)
  ε > 0 || return 0.0

  _Etot = Etot === nothing ? Float64(E0) + ε : Float64(Etot)

  # -- find s_ε
  roots_ε = find_s_all_Eres(pathfit, ε; tol=tol)
  isempty(roots_ε) && return 0.0
  s_ε = select_root(roots_ε; branch=branch_ε)
  s_ε === nothing && return 0.0

  # -- find S_E
  roots_Etot = find_s_all_U(pathfit, _Etot; tol=tol)
  isempty(roots_Etot) && return 0.0

  s_Etot = if branch_Etot === :same_side
    select_root(roots_Etot; branch=_branch_from_s(s_ε))
  elseif branch_Etot === :closest_to_sε
    select_root(roots_Etot; branch=:closest, sref=s_ε)
  else
    select_root(roots_Etot; branch=branch_Etot)
  end
  s_Etot === nothing && return 0.0

  Γε = Γ(pathfit, s_ε)
  slope = abs(dUds(pathfit, s_Etot))
  slope == 0.0 && return 0.0

  ζ2 = abs2(ζ_eval(s_Etot)) # |ζ(s_Etot)|²
  k2 = 2.0 * ε # k² = 2E in au

  return 2π^2 * Γε * ζ2 / (k2 * slope)

end

"""
    σ(pathfit, energies; Etot = nothing, tol=1e-10)

Compute the dissociation cross section for an array of electron energies `energies`
"""
function σ(
      pathfit :: DissociationPathwayFit
    , energies :: AbstractVector{<:Real}
    ; kwargs...
  )
  return [σ(pathfit, ε; kwargs...) for ε in energies]
end

function compute_partial_dissociation(
      pathfit :: DissociationPathwayFit
    ; channel :: Channel
    , Egrid :: AbstractVector{<:Real}
    , mode_freqs :: AbstractDict{String, <:Real}
    , targets_by_mode :: AbstractDict{String, <:TargetState}
    , spinmult :: Union{Int, Nothing} = nothing
    , target_spinmult :: Union{Int, Nothing} = nothing
    , spin_mode :: Symbol = :auto
    , kwargs...
  )

  E0 = channel_E0(channel, mode_freqs)
  ζ_eval = ζ_path_fit(pathfit, targets_by_mode; kind=:linear)

  family = _family_name(pathfit.channel)
  spinmult, spinwght = resolve_spinwght(family, spinmult, target_spinmult, spin_mode)

  σvals = σ(pathfit, Egrid; E0=E0, ζ_eval=ζ_eval, kwargs...)
  σvals .*= spinwght

  return PartialDissociationResult(
      pathfit.channel
    , family
    , spinmult
    , spinwght
    , copy(Float64.(Egrid))
    , σvals
    , pathfit.path
  )

end

"""
    compute_total_dissociation(pathfits, spec; tol=1e-10)

Compute all partial dissociation cross sections and combine them to form the total cross sections,
and determine their relative contributions.
"""
function compute_total_dissociation(
      pathfits :: AbstractVector{<:DissociationPathwayFit}
    ; channels :: AbstractVector{Channel}
    , Egrid :: AbstractVector{<:Real}
    , mode_freqs :: AbstractDict{String, <:Real}
    , targets_by_mode :: AbstractDict{String, <:TargetState}
    , channel_spinmults :: AbstractDict{String, <:Integer} = Dict{String, Int}()
    , target_spinmult :: Union{Int, Nothing} = nothing
    , spin_mode :: Symbol = :auto
    , kwargs...
  )

  nchannels = length(channels)
  npathfits = length(pathfits)
  nchannels == npathfits || error("Number of channels ($nchannels) does not match number of pathfits ($npathfits) ")

  partials = [
    compute_partial_dissociation(
        pathfits[i]
      ; channel=channels[i]
      , Egrid=Egrid
      , mode_freqs=mode_freqs
      , targets_by_mode=targets_by_mode
      , spinmult=get(channel_spinmults, pathfits[i].channel, nothing)
      , target_spinmult=target_spinmult
      , spin_mode=spin_mode
      , kwargs...
    )
    for i in eachindex(pathfits, channels)
  ]

  isempty(partials) && error("No partial results to combine !")

  energies = copy(Float64.(Egrid))
  σtot = zeros(Float64, length(energies))

  for part in partials
    σtot .+= part.σ
  end

  channel_fractions = Dict{String, Vector{Float64}}()

  for part in partials
    frac = similar(part.σ)
    for i in eachindex(frac)
      frac[i] = σtot[i] > 0 ? part.σ[i] / σtot[i] : 0.0
    end
    channel_fractions[part.channel] = frac
  end

  familyσ = Dict{String, Vector{Float64}}()
  for part in partials
    if !haskey(familyσ, part.family)
      familyσ[part.family] = zeros(Float64, length(energies))
    end
    familyσ[part.family] .+= part.σ
  end

  family_fractions = Dict{String, Vector{Float64}}()
  for (fam, σfam) in familyσ
    frac = similar(σfam)
    for i in eachindex(frac)
      frac[i] = σtot[i] > 0 ? σfam[i] / σtot[i] : 0.0
    end
    family_fractions[fam] = frac
  end

  return TotalDissociationResult(
      partials
    , energies
    , σtot
    , channel_fractions
    , family_fractions
  )

end
