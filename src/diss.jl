using Statistics: mean
using LinearAlgebra: norm

###############
### HELPERS ###
###############

"""
    _require_col(table, name)

Return the column index for `name` or throw an error
"""
function _require_col(table :: Table, name :: AbstractString)
  haskey(table.colmap, String(name)) || error("Table is missing the required column '$name'")
  return table.colmap[String(name)]
end

"""
    _require_modefreq(spec, mode_name)

Return the mode frequency for `mode_name` or throw an error
"""
function _require_modefreq(spec :: DissCalcSpec, mode_name :: AbstractString)
  haskey(spec.mode_freqs, String(mode_name)) || error("No mode frequency provided for mode '$mode_name'")
  return Float64(spec.mode_freqs[String(mode_name)])
end

"""
    _channel_mode_pairs(spec)

Flattened channel and mode iterator -> (Channel Name, Mode)
"""
function _channel_mode_pairs(spec :: DissCalcSpec)
  pairs = Tuple{String, ChannelMode}[]
  for ch in spec.channels
    for cm in ch.modes
      push!(pairs, (ch.name, cm))
    end
  end
  return pairs
end

"""
    _surface_modes(channel)

Mode names in coordinate order for the channel surface
"""
_surface_modes(channel :: Channel) = [cm.mode_name for cm in channel.modes]

"Reference value of Eres at Q=0"
_fit_Eres0(fit :: ModeFit) = fit.Eres_fit(0.0)
"Reference value of Γ at Q=0"
_fit_Γ0(fit :: ModeFit) = fit.Γ_fit(0.0)

"Converts `q` to a Float64 and makes sure that it has the right length"
function _qvec(surface :: ChannelSurface, q :: AbstractVector)
  length(q) == length(surface.modes) || error(
    "Vector q has length $(length(q)), expected $(length(surface.modes))"
  )
  return Float64.(q)
end

"Simple finite difference [f(x+h) - f(x-h)] / 2h"
function _fd(f, x; h :: Float64 = 1e-4)
  xp = Float64(x) + h
  xm = Float64(x) - h
  return (f(xp) - f(xm)) / 2h
end

"Derivative wrapper for _fd(f, x;h=1e-4)"
_deriv(f, x :: Real; h :: Float64 = 1e-4) = _fd(f, x; h=h)

"""
    _deriv(f, x)

Evaluates the analytic derivative of a polynomial function
"""
function _deriv(f :: PolyInterp1D, x :: Real)
  x = Float64(x)
  f.clampx && (x = clamp(x, f.xmin, f.xmax))
  n = length(f.coeffs)
  n <=1 && return 0.0 # derivative of a constant
  # -- horner derivative
  y = (n-1)*f.coeffs[end]
  for k in (n-1):-1:2
    y = muladd(y, x, (k-1)*f.coeffs[k])
  end
  return y
end

###################
###################
### P U B L I C ###
##VVVVVVVVVVVVVVV##
#vvvvvvvvvvvvvvvvv#

"""
    prepare_dissociation_model(widetable, spec; targets_by_mode, kwargs...)

Widetable -> dissociation model pipeline
"""
function prepare_dissociation_model(
    widetable :: Table
  , spec :: DissCalcSpec
  ; targets_by_mode :: AbstractDict{String, <:Any} = Dict{String, Any}()
  , kwargs...
  )
  validate_calc_spec(widetable, spec)
  longtable = widetable_to_longtable(widetable, spec)
  curves = parse_mode_curves(longtable; qcol = spec.qcol)
  fits = fit_mode_curves(curves, spec; targets_by_mode = targets_by_mode)
  surfaces = build_channel_surfaces(fits, spec)
  paths = build_dissociation_paths(surfaces, spec; kwargs...)
  return (
      longtable = longtable
    , curves = curves
    , fits = fits
    , surfaces = surfaces
    , paths = paths
  )
end

###############
### PARSING ###
###############

"""
    validate_calc_spec(widetable, spec)

Given a widetable with resonance positions and widths for a set of channels in their
respective modes, validate that the table is compatible with `DissCalcSpec` at all
and throw errors otherwise.

Checks the following :
- `qcol` exists
- channel names are unique
- mode names are unique for each channel
- required `'label'_E` and `'label'_Γ` columns exist
- if `V_kind === :harmonic`, then mode frequencies must also exist
"""
function validate_calc_spec(widetable :: Table, spec :: DissCalcSpec)

  _require_col(widetable, spec.qcol)

  seen_channels = Set{String}()

  for ch in spec.channels
    ch.name in seen_channels && error("Duplicate channel name '$(ch.name)' in spec !")
    push!(seen_channels, ch.name)

    isempty(ch.modes) && error("Channel '$(ch.name)' has no modes !")

    seen_modes = Set{String}()
    for cm in ch.modes
      cm.mode_name in seen_modes && error("Channel '$(ch.name)' repeats mode '$(cm.mode_name)' !")
      push!(seen_modes, cm.mode_name)

      _require_col(widetable, "$(cm.label)_E")
      _require_col(widetable, "$(cm.label)_Γ")
      spec.interp.V_kind === :harmonic && _require_modefreq(spec, cm.mode_name)

    end

  end

  return spec

end

"""
    parse_mode_curves(longtable; qcol="Q")


Parses a long calculation table with columns

    Q, channel, mode, label, Eres, Γ

into a vector of `ModeCurve`s, skipping rows with missing `Q`, `Eres`, or `Γ`.
"""
function parse_mode_curves(longtable :: Table; qcol :: String = "Q")

  Q_idx    = _require_col(longtable, qcol)
  ch_idx   = _require_col(longtable, "channel")
  mode_idx = _require_col(longtable, "mode")
  lbl_idx  = _require_col(longtable, "label")
  E_idx    = _require_col(longtable, "Eres")
  Γ_idx    = _require_col(longtable, "Γ")

  groups = Dict{
      Tuple{String,          String,          String}
    , Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
   }()

  for i in axes(longtable.rows, 1)
    Q    = longtable.rows[i, Q_idx]
    ch   = longtable.rows[i, ch_idx]
    mode = longtable.rows[i, mode_idx]
    lbl  = longtable.rows[i, lbl_idx]
    E    = longtable.rows[i, E_idx]
    Γ    = longtable.rows[i, Γ_idx]

    (ismissing(Q) || ismissing(E) || ismissing(Γ)) && continue

    key = (String(ch), String(mode), String(lbl))
    haskey(groups, key) || (groups[key] = (Float64[], Float64[], Float64[]))

    Qv, Ev, Γv = groups[key]

    push!(Qv, Float64(Q))
    push!(Ev, Float64(E))
    push!(Γv, Float64(Γ))
  end

  curves = ModeCurve[]
  for ((ch, mode, lbl), (Qv, Ev, Γv)) in groups
    idx = sortperm(Qv)
    push!(curves, ModeCurve(
        ch
      , mode
      , lbl
      , Qv[idx]
      , Ev[idx]
      , Γv[idx]
    ))
  end

  sort!(curves, by = c -> (c.channel_name, c.mode_name, c.label))

  return curves

end

"""
    parse_mode_curves(widetable, spec)

Convenience wrapper for directly passing a widetable

    widetable -> longtable -> ModeCurves
"""
function parse_mode_curves(widetable :: Table, spec :: DissCalcSpec)
  longtable = widetable_to_longtable(widetable, spec)
  return parse_mode_curves(longtable; qcol = spec.qcol)
end

###############
### FITTING ###
###############

"""
    fit_mode_curve

Builds a `ModeFit` from a `ModeCurve`, constructing callable fits for the following :
- `Eres(Q)`
- `Γ(Q)`
- `V(Q)`

Arguments
---------
- `curve` :
    Resonance curve extracted from the longtable
- `interp` :
    `InterpSpec` describing how to fit/interpolate `Eres`, `Γ`, and `V`
- `target` :
    Optional `TargetState` for this mode; required when `interp.V_kind` needs the target potential curve
- `mode_freqs` :
    Optional mode frequencies; needed when `interp.V_kind` requries them,
    e.g., `interp_Vkind === :harmonic`
"""
function fit_mode_curve(
    curve :: ModeCurve
  ; interp :: InterpSpec
  , target :: Union{TargetState, Nothing} = nothing
  , mode_freqs :: AbstractDict{String, <:Real} = Dict{String, Float64}()
  )

  isempty(curve.Q) &&
    error("Can't fit empty ModeCurve for channel '$(curve.channel_name)', 'mode = $(curve.mode_name)'")

  Qmin, Qmax = extrema(curve.Q)

  # -- Eres(Q)
  Eres_fit = _make_interp(curve.Q, curve.Eres; kind=interp.Eres_kind)

  # -- Γ(Q)
  Γ_fit = _make_interp(curve.Q, curve.Γ; kind=interp.Γ_kind)

  # -- V(Q)
  V_fit = if interp.V_kind === :harmonic

    haskey(mode_freqs, curve.mode_name) || error("Need mode frequency for mode '$(curve.mode_name)' when V_kind = :harmonic")
    ω = Float64(mode_freqs[curve.mode_name])

    Qeq = if target === nothing
      0.0
    else
      target.Q[argmin(target.V)]
    end

    # -- dimensionless coordinate ! V(Q) = ωQ²/2, NOT ω²Q²/2
    Vref0 = 0.5 * ω * (0.0 - Qeq)^2
    Q -> 0.5 * ω * (Float64(Q) - Qeq)^2 - Vref0

  elseif interp.V_kind === :linear

    target === nothing && error("Need a TargetState for mode '$(curve.mode_name)' when V_kind = :linear")

    # -- linear interpolation s.t. V(0) = 0.0
    V0 = linterp(target.Q, target.V, 0.0)
    _make_interp(target.Q, target.V .- V0; kind=:linear)

  elseif interp.V_kind === :poly2

    target === nothing && error("Need a TargetState for mode '$(curve.mode_name)' when V_kind = :poly2")

    # -- quadratic fit of the target potential curve s.t. V(0)=0
    V0 = linterp(target.Q, target.V, 0.0)
    _make_interp(target.Q, target.V .- V0; kind=:poly2)

  else
    error("Unsupported V_kind '$(interp.V_kind)'")
  end

  return ModeFit(
      curve
    , target
    , Eres_fit
    , Γ_fit
    , V_fit
    , Qmin
    , Qmax
  )

end

"""
    fit_mode_curves(curves, spec; targets_by_mode)

Fit a whole vector of `ModeCurve`s into a vector of `ModeFit`s.
`targets_by_mode` should map mode_name => TargetState (because V(Q) in this sense
is a property of the normal mode, not of the channel)
"""
function fit_mode_curves(
    curves :: AbstractVector{ModeCurve}
  , spec :: DissCalcSpec
  ; targets_by_mode :: AbstractDict{String, <:Any} = Dict{String, Any}()
  )

  fits = ModeFit[]

  for curve in curves
    target = get(targets_by_mode, curve.mode_name, nothing)

    push!(fits, fit_mode_curve(
        curve
      ; interp = spec.interp
      , target = target
      , mode_freqs = spec.mode_freqs
    ))
  end

  return fits

end

"""
   fits_for_channel(fits, channel)

Extracts the `ModeFit`s belonging to the channel `channel`, keyed by mode name.
Matches on channel name, mode name, and label
"""
function fits_for_channel(
    fits :: AbstractVector{ModeFit}
  , channel :: Channel
  )

  out = Dict{String, ModeFit}()

  for cm in channel.modes
    idx = findfirst(f ->
      f.curve.channel_name == channel.name &&
      f.curve.mode_name == cm.mode_name &&
      f.curve.label == cm.label
    , fits)

    idx === nothing && error("No ModeFit found for channel='$(channel.name)', mode='$(cm.mode_name)', label='$(cm.label)'")

    out[cm.mode_name] = fits[idx]

  end

  return out

end

#############################
### BUILDING THE SURFACES ###
#############################

"""
    build_channel_surface(fits, channel, spec)

Buils a `ChannelSurface` for one channel from the fitted mode curves
"""
function build_channel_surface(
    fits :: AbstractVector{ModeFit}
  , channel :: Channel
  , spec :: DissCalcSpec
  )

  fitdict = fits_for_channel(fits, channel)
  modes   = _surface_modes(channel)

  isempty(modes) && error("Channel '$(channel.name)' has no modes !")

  # -- reference values at Q=0, averaged over mode fits
  E0arr = Float64[_fit_Eres0(fitdict[m]) for m in modes]
  Γ0arr = Float64[_fit_Γ0(fitdict[m])    for m in modes]
  Eres0 = mean(E0arr)
  Γ0 = mean(Γ0arr)

  return ChannelSurface(
      channel.name
    , modes
    , fitdict
    , spec.mode_freqs
    , spec.combine_Eres
    , spec.combine_Γ
    , Eres0
    , Γ0
  )

end

"""
    build_channel_surfaces(fits, spec)

Build one `ChannelSurface` per channel in the calculation spec
"""
function build_channel_surfaces(
    fits :: AbstractVector{ModeFit}
  , spec :: DissCalcSpec
  )
  surfaces = ChannelSurface[]
  for ch in spec.channels
    push!(surfaces, build_channel_surface(fits, ch, spec))
  end
  return surfaces
end

###############################
### EVALUATING THE SURFACES ###
###############################

"""
    Eres(surface, q)

Evaluates the channel resonance energy at the coordinate vector q for each mode simultaneously:
    Eres(q) = Eres0 + Σ_i [Eres_i(q_i) - Eres_i(0)]
"""
function Eres(
      surface :: ChannelSurface
    , q :: AbstractVector
  )

  qv = _qvec(surface, q)

  surface.combine_Eres !== :additive && error("Unsupported combine_Eres ($(surface.combine_Eres))")

  val = surface.Eres0
  for (i, mode) in enumerate(surface.modes)
    fit = surface.fits[mode]
    val += fit.Eres_fit(qv[i]) - _fit_Eres0(fit)
  end
  return val

end

"""
    Γ(surface, q)

Evaluate the channel resonance width at the coordinate vector q for each mode simultaneously:
- surface.combine_Γ = :additive -> Γ(q) = Γ0 + Σ_i [Γ_i(q_i) - Γ_i(0)]
- surface.combine_Γ = :max -> Γ(q) = maximum(Γ_i(q_i))
"""
function Γ(
      surface :: ChannelSurface
    , q :: AbstractVector
  )

  qv = _qvec(surface, q)

  if surface.combine_Γ === :additive

    val = surface.Γ0
    for (i, mode) in enumerate(surface.modes)
      fit = surface.fits[mode]
      val += fit.Γ_fit(qv[i]) - _fit_Γ0(fit)
    end
    return max(0.0, val)

  elseif surface.combine_Γ === :max

    vals = Float64[]
    for (i, mode) in enumerate(surface.modes)
      fit = surface.fits[mode]
      push!(vals, fit.Γ_fit(qv[i]))
    end
    return max(0.0, maximum(vals))

  else
    error("Unsupported combine_Γ = '$(surface.combine_Γ)'")
  end

end

"""
    V(surface, q)

Evaluates the target molecule's potential energy at coordinate vector q.
Assumes each `ModeFit.V_fit` is already referenced s.t. V_i(0) = 0:
    V(q) = Σ_i V_i(q_i)
"""
function V(
    surface :: ChannelSurface
  , q :: AbstractVector
  )

  qv = _qvec(surface, q)

  val = 0.0
  for (i, mode) in enumerate(surface.modes)
    fit = surface.fits[mode]
    val += fit.V_fit(qv[i])
  end

  return val

end

#################
### GRADIENTS ###
#################

"""
    U(surface, q)

Evaluates the dissociative surface U(q) = V(q) + Eres(q)
"""
U(surface :: ChannelSurface, q :: AbstractVector) = V(surface, q) + Eres(surface, q)


"""
    ∇E(surface, q)

Computes the gradient of Eres with respect to channel coordinates
"""
function ∇Eres(
    surface :: ChannelSurface
  , q :: AbstractVector
  )

  qv = _qvec(surface, q)
  ∇E = zeros(Float64, length(qv))

  for (i, mode) in enumerate(surface.modes)
    fit = surface.fits[mode]
    ∇E[i] = _deriv(fit.Eres_fit, qv[i])
  end

  return ∇E

end

"""
    ∇V(surface, q)

Computes the gradient of V with respect to channel coordinate
"""
function ∇V(
    surface :: ChannelSurface
  , q :: AbstractVector
  )

  qv = _qvec(surface, q)
  ∇V = zeros(Float64, length(qv))

  for (i, mode) in enumerate(surface.modes)
    fit = surface.fits[mode]
    ∇V[i] = _deriv(fit.V_fit, qv[i])
  end

  return ∇V

end

"""
    ∇U(surface, q)

Computes the gradient of U = V + Eres with respect to channel coordinate
"""
∇U(surface :: ChannelSurface, q :: AbstractVector) = ∇V(surface, q) + ∇Eres(surface, q)

"""
   dUds(surface, q)

Calculates the derivative of U with respect to the steepest descent coordinate s:
   dq/ds = -∇U/|∇U| (vector quantity because q is a vector)
   dU/ds = ∇U · dq/ds = -|∇U|
"""
dUds(surface :: ChannelSurface, q :: AbstractVector) = -norm(∇U(surface, q))

############################
### DISSOCIATION PATHWAY ###
############################

"""
    is_in_fit_range(surface, q; pad=0.0)

Checks whether the coodinates q are within the surface boundaries, with some optional padding
"""
function is_in_fit_range(
    surface :: ChannelSurface
  , q :: AbstractVector
  ; pad :: Float64 = 0.0
  )
  qv = _qvec(surface, q)
  for (i, mode) in enumerate(surface.modes)
    fit = surface.fits[mode]
    (qv[i] < fit.qmin - pad || qv[i] > fit.qmax + pad) && return false
  end
  return true
end

function build_dissociation_path(
    surface :: ChannelSurface
  , spec :: DissCalcSpec
  ; q0 :: AbstractVector = zeros(length(surface.modes))
  , range_pad :: Float64 = 0.0
  , grad_tol :: Float64 = 1e-8
  )

  q = _qvec(surface, q0)
  length(q) == length(surface.modes) || error("q0 length ($(length(q0))) must match `surface.modes` length")

  sval = 0.0

  svec = Float64[]
  qrows = Vector{Vector{Float64}}()
  Evec = Float64[]
  Γvec = Float64[]
  Vvec = Float64[]
  Uvec = Float64[]
  dUdsvec = Float64[]

  while true
    is_in_fit_range(surface, q; pad=range_pad) || break

    g = ∇U(surface, q)
    gnorm = norm(g)
    gnorm <= grad_tol && break

    E = Eres(surface, q)
    _Γ = Γ(surface, q)
    Vq = V(surface, q)
    Uq = Vq + E
    _dUds = -gnorm

    push!(svec, sval)
    push!(qrows, copy(q))
    push!(Evec, E)
    push!(Γvec, _Γ)
    push!(Vvec, Vq)
    push!(Uvec, Uq)
    push!(dUdsvec, _dUds)

    sval >= spec.path.smax && break

    q .-= spec.path.ds .* g ./ gnorm
    sval += spec.path.ds

  end

  qmat = Matrix{Float64}(undef, length(qrows), length(surface.modes))
  for i in eachindex(qrows)
    qmat[i, :] .= qrows[i]
  end

  return DissociationPathway(
    surface.channel_name
    , copy(surface.modes)
    , svec
    , qmat
    , Evec
    , Γvec
    , Vvec
    , Uvec
    , dUdsvec
  )

end

"""
    build_dissociation_paths(surface, spec; kwargs...)

Wrapper for multiple dissociation pathways.
"""
function build_dissociation_paths(
    surfaces :: AbstractVector{ChannelSurface}
  , spec :: DissCalcSpec
  ; kwargs...
  )
  pathways = DissociationPathway[]
  for surface in surfaces
    push!(pathways, build_dissociation_path(surface, spec; kwargs...))
  end
  return pathways
end
