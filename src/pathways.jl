using Statistics: mean
using LinearAlgebra: norm

export prepare_dissociation_model
export validate_calc_spec
export parse_mode_curves
export fit_mode_curve, fit_mode_curves
export fits_for_channel
export build_channel_surface, build_channel_surfaces
export fit_path, fit_paths
export Eres, Γ, V, U, dUds
export build_dissociation_path,  build_dissociation_path_autostart
export build_dissociation_paths, build_dissociation_paths_autostart

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
    _require_s_in_range(pathfit, s)

Return the dissociation coordinate s or throw an error if it is not in the expected bounds
"""
function _require_s_in_range(pathfit :: DissociationPathwayFit, s :: Real)
  smin, smax = extrema(pathfit.s)
  sval = Float64(s)
  smin <= sval <= smax || error("Requested s = $sval is outside the fitted range [$smin,$smax]")
  return sval
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
_fit_Eres0(fit :: ModeFit) = fit.Eres_fit(clamp(0.0, fit.qmin, fit.qmax))
"Reference value of Γ at Q=0"
_fit_Γ0(fit :: ModeFit) = fit.Γ_fit(clamp(0.0, fit.qmin, fit.qmax))

"Converts `q` to a Float64 and makes sure that it has the right length"
function _qvec(surface :: ChannelSurface, q :: AbstractVector)
  length(q) == length(surface.modes) || error(
    "Vector q has length $(length(q)), expected $(length(surface.modes))"
  )
  return Float64.(q)
end

"Simple finite difference [f(x+h) - f(x-h)] / 2h"
@inline function _fd(f, x :: Float64; h :: Float64 = 1e-4)

  x0 = Float64(x)
  xp = x0 + h
  xm = x0 - h

  y0 = f(x0)

  okp, okm = (true, true)
  yp,  ym  = (0.0,  0.0)

  # -- get f(x) if in bounds
  try
    yp = f(xp)
  catch
    okp = false
  end

  # -- get f(x) if in bounds
  try
    ym = f(xm)
  catch
    okm = false

  end

  if okp && okm
    return (yp - ym) / (2h)
  elseif okp
    return (yp - y0) / h
  elseif okm
    return (y0 - ym) / h
  else
    error("Could not evaluate finite-difference derivative at x=$x0; both evaluations failed for
          xm = $xm
          xp = $xp")
  end
end
@inline _fd(f, x::Real; h::Float64 = 1e-4) = _fd(f, Float64(x); h=h)

"Derivative wrapper for _fd(f, x;h=1e-4)"
@inline _deriv(f, x :: Real;    h :: Float64 = 1e-4) = _fd(f, Float64(x); h=h)
@inline _deriv(f, x :: Float64; h :: Float64 = 1e-4) = _fd(f, x; h=h)


"""
    _deriv(f, x)

Evaluates the analytic derivative of a polynomial function
"""
@inline function _deriv(f :: PolyInterp1D, x :: Float64)
  xx = f.clampx ? clamp(x, f.xmin, f.xmax) : x
  n = length(f.coeffs)
  n <=1 && return 0.0 # derivative of a constant
  # -- horner derivative
  y = (n-1)*f.coeffs[end]
  @inbounds for k in (n-1):-1:2
    y = muladd(y, xx, (k-1)*f.coeffs[k])
  end
  return y
end
@inline _deriv(f :: PolyInterp1D, x::Real) = _deriv(f, Float64(x))

@inline function _deriv(f :: LinearInterp1D, x::Float64)
  xs, ys = f.x, f.y
  n = length(xs)
  n >=2 || return 0.0
  if x<= xs[1]
    i=1
  elseif x>= xs[end]
    i=n-1
  else
    i=searchsortedlast(xs,x)
    i=clamp(i, 1, n-1)
  end
  dx, dy = xs[i+1] - xs[i], ys[i+1] - ys[i]
  return dy/dx
end
@inline _deriv(f :: LinearInterp1D, x :: Real) = _deriv(f, Float64(x))

@inline function _deriv(h :: HarmonicPotential, Q :: Float64)
  return h.ω * (Q - h.Qeq)
end
@inline _deriv(h :: HarmonicPotential, Q :: Real) = _deriv(h, Float64(Q))

@inline function _deriv(f :: Union{PCHIPInterp1D, AkimaInterp1D}, x :: Float64; h :: Float64 = 1e-4)
  xx = f.clampx ? clamp(x, f.xmin, f.xmax) : x
  xl = max(f.xmin, xx - h)
  xr = min(f.xmax, xx + h)
  xr == xl && return 0.0
  return (f(xr) - f(xl)) / (xr - xl)
end
@inline _deriv(f :: Union{PCHIPInterp1D, AkimaInterp1D}, x :: Real; h :: Float64 = 1e-4) = _deriv(f, Float64(x); h=h)

"""
    _find_brackets(x, y, ytarget)

Finds intervals `[x[i], x[i+1]]` where `(y - ytarget)` changes sign.
Returns a vector of index pairs `(i,i+1)`
"""
function _find_brackets(
    x :: AbstractVector{<:Real}
  , y :: AbstractVector{<:Real}
  , ytarget :: Real
  )

  n = length(x)
  n == length(y) || error("x/y length mismatch")
  n < 2 && return Tuple{Int, Int}[]

  out = Tuple{Int, Int}[]

  for i in 1:(n-1)
    f1 = y[i]   - ytarget
    f2 = y[i+1] - ytarget
    if f1 == 0.0
      push!(out, (i,i))
    elseif f1*f2 < 0.0 || f2 == 0.0
      push!(out, (i,i+1))
    end
  end

  return out

end

"""
    _bisect_root(f, a, b; tol=1e=10, maxiter=100)

Simple bisection root finder for a function f in the interval [a,b]
"""
function _bisect_root(
    f
  , a :: Real
  , b :: Real
  ; tol :: Float64 = 1e-10
  , maxiter :: Int = 100
  )

  fa, fb = f.((a,b))

  fa == 0.0 && return Float64(a)
  fb == 0.0 && return Float64(b)

  fa * fb > 0.0 && error("Root is not in [$a, $b]")

  lo, hi   = Float64.((a, b))
  flo, fhi = fa, fb

  for _ in 1:maxiter
    mid = 0.5 * (lo+hi)
    fm = f(mid)

    abs(fm)      <= tol && return mid
    abs(hi - lo) <= tol && return mid

    if flo*fm <= 0.0
      hi = mid
      fhi = fm
    else
      lo = mid
      flo = fm
    end
  end

  return 0.5 * (lo+hi)

end

function _merge_partial_dissociation_paths(
    neg :: DissociationPathway
  , pos :: DissociationPathway
  )

  neg.channel == pos.channel || error("Channel mismatch in path merge")
  neg.modes   == pos.modes   || error("Mode mismatch in path merge")

  # -- reverse negative branch so that s is ascending
  s_neg    = reverse(neg.s)
  q_neg    = reverse(neg.q; dims=1)
  E_neg    = reverse(neg.Eres)
  Γ_neg    = reverse(neg.Γ)
  V_neg    = reverse(neg.V)
  U_neg    = reverse(neg.U)
  dUds_neg = reverse(neg.dUds)

  # -- remove duplicated s=0 from positive branch if it's in both
  istart = (!isempty(s_neg) && !isempty(pos.s) && s_neg[end] == 0.0 && pos.s[1] == 0.0) ? 2 : 1

  s_all    = vcat(s_neg,    pos.s[istart:end])
  q_all    = vcat(q_neg,    pos.q[istart:end, :])
  E_all    = vcat(E_neg,    pos.Eres[istart:end])
  Γ_all    = vcat(Γ_neg,    pos.Γ[istart:end])
  V_all    = vcat(V_neg,    pos.V[istart:end])
  U_all    = vcat(U_neg,    pos.U[istart:end])
  dUds_all = vcat(dUds_neg, pos.dUds[istart:end])

  return DissociationPathway(
      pos.channel
    , copy(pos.modes)
    , s_all
    , q_all
    , E_all
    , Γ_all
    , V_all
    , U_all
    , dUds_all
  )
end

"Returns a small list of candidate starting points `q0` for dissociation path construction"
function _get_default_q0_candidates(surface :: ChannelSurface)

  n = length(surface.modes)

  # -- try the original default
  cands = Vector{Vector{Float64}}()
  push!(cands, zeros(n))

  # -- try nearest in range for each mode fit
  qclamp0 = [clamp(0.0, fit.qmin, fit.qmax) for fit in surface.fits]
  qclamp0 == cands[1] || push!(cands, qclamp0)

  # -- try lower/upper corners of fit box
  qmins = [fit.qmin for fit in surface.fits]
  qmaxs = [fit.qmax for fit in surface.fits]
  push!(cands, qmins)
  push!(cands, qmaxs)

  # -- try modewise midpoints
  qmid = [(fit.qmin + fit.qmax)/2 for fit in surface.fits]
  push!(cands, qmid)

  return cands

end

# @inline _surface_interp_override(
#       surface_interp_overrides :: AbstractDict{String, <:Any}
#     , channel_name :: AbstractString
#    ) = get(surface_interp_overrides, String(channel_name), nothing)

# @inline _surface_extrap_override(
#       surface_extrap_overrides :: AbstractDict{String, <:Any}
#     , channel_name :: AbstractString
#    ) = get(surface_extrap_overrides, String(channel_name), nothing)

@inline function _merge_res_interp_spec(
      interp :: InterpSpec
    , override
  )
  override == nothing && return interp
  return InterpSpec(
      isnothing(override.Eres_kind) ? interp.Eres_kind : override.Eres_kind
    , isnothing(override.Γ_kind)    ? interp.Γ_kind    : override.Γ_kind
    , interp.V_kind
    , interp.path_kind
  )
end

@inline function _merge_res_extrap_spec(
      extrap :: ExtrapSpec
    , override
  )
  override == nothing && return extrap
  return ExtrapSpec(
      isnothing(override.Eres_left_kind)   ? extrap.Eres_left_kind   : override.Eres_left_kind
    , isnothing(override.Eres_right_kind)  ? extrap.Eres_right_kind  : override.Eres_right_kind
    , isnothing(override.Γ_left_kind)      ? extrap.Γ_left_kind      : override.Γ_left_kind
    , isnothing(override.Γ_right_kind)     ? extrap.Γ_right_kind     : override.Γ_right_kind
    , isnothing(override.EresΓ_Qpad_left)  ? extrap.EresΓ_Qpad_left  : override.EresΓ_Qpad_left
    , isnothing(override.EresΓ_Qpad_right) ? extrap.EresΓ_Qpad_right : override.EresΓ_Qpad_right
    , extrap.V_left_kind
    , extrap.V_right_kind
    , extrap.V_Qpad_left
    , extrap.V_Qpad_right
  )
end

# @inline function _merge_interp_spec(
#       interp :: InterpSpec
#     , override
#   )

#   override === nothing && return interp

#   return InterpSpec(
#       isnothing(override.Eres_kind) ? interp.Eres_kind : override.Eres_kind
#     , isnothing(override.Γ_kind)    ? interp.Γ_kind    : override.Γ_kind
#     , isnothing(override.V_kind)    ? interp.V_kind    : override.V_kind
#     , isnothing(override.path_kind) ? interp.path_kind : override.path_kind
#   )
# end

# @inline function _merge_extrap_spec(
#       extrap :: ExtrapSpec
#     , override
#   )

#   override === nothing && return extrap

#   return ExtrapSpec(
#       isnothing(override.Eres_left_kind)   ? extrap.Eres_left_kind   : override.Eres_left_kind
#     , isnothing(override.Eres_right_kind)  ? extrap.Eres_right_kind  : override.Eres_right_kind
#     , isnothing(override.Γ_left_kind)      ? extrap.Γ_left_kind      : override.Γ_left_kind
#     , isnothing(override.Γ_right_kind)     ? extrap.Γ_right_kind     : override.Γ_right_kind
#     , isnothing(override.EresΓ_Qpad_left)  ? extrap.EresΓ_Qpad_left  : override.EresΓ_Qpad_left
#     , isnothing(override.EresΓ_Qpad_right) ? extrap.EresΓ_Qpad_right : override.EresΓ_Qpad_right
#     , isnothing(override.V_left_kind)      ? extrap.V_left_kind      : override.V_left_kind
#     , isnothing(override.V_right_kind)     ? extrap.V_right_kind     : override.V_right_kind
#     , isnothing(override.V_Qpad_left)      ? extrap.V_Qpad_left      : override.V_Qpad_left
#     , isnothing(override.V_Qpad_right)     ? extrap.V_Qpad_right     : override.V_Qpad_right
#   )
# end

# @inline function _path_kind_with_override(
#       default_kind :: Symbol
#     , override
#   )
#   override === nothing && return default_kind
#   return isnothing(override.path_kind) ? default_kind : override.path_kind
# end

"""
    _edge_poly2_interp(x, y; side=:left, nedge=5)

Quadratic fit to edge points, used for local edge extrapolation.
Falls back to `nothing` if <3 points are available.
"""
function _edge_poly2_interp(
      x :: AbstractVector{<:Real}
    , y :: AbstractVector{<:Real}
    ; side=:left
    , nedge=5
  )

  xs, ys = _sorted_xy(x, y)
  n = length(xs)
  n < 3 && return nothing

  nuse = min(nedge, n)
  if side === :left
    idx = 1:nuse
  elseif side === :right
    idx = (n - nuse + 1):n
  else
    throw(ArgumentError("side must be :left or :right"))
  end

  length(idx) < 3 && return nothing
  c = _polyfit(xs[idx], ys[idx], 2)
  return PolyInterp1D(c, xs[first(idx)], xs[last(idx)], false)

end

function _wrap_with_extrapolation(
      f_in
    , xs :: AbstractVector{<:Real}
    , ys :: AbstractVector{<:Real}
    ; left_kind :: Symbol = :none
    , right_kind :: Symbol = :none
    , Qpad_left :: Real = 0.0
    , Qpad_right :: Real = 0.0
    , nedge_left :: Int = 5
    , nedge_right :: Int = 5
    , harmonic_model = nothing
    , clamp_nonnegative :: Bool = false
  )

  xfit, yfit = _sorted_xy(xs, ys)
  xmin, xmax = extrema(xfit)

  yL, yR = f_in.((xmin, xmax))

  mL = _deriv(f_in, xmin)
  mR = _deriv(f_in, xmax)

  # -- edge secants from raw data
  sL = length(xfit) >= 2 ? (yfit[2]   - yfit[1])     / (xfit[2]   - xfit[1])     : mL
  sR = length(xfit) >= 2 ? (yfit[end] - yfit[end-1]) / (xfit[end] - xfit[end-1]) : mR

  polyL = left_kind  === :poly2_edge ? _edge_poly2_interp(xfit, yfit; side=:left,  nedge=nedge_left)  : nothing
  polyR = right_kind === :poly2_edge ? _edge_poly2_interp(xfit, yfit; side=:right, nedge=nedge_right) : nothing

  @inline _maybe_nonneg(y) = clamp_nonnegative ? max(0.0, y) : y

  function fout(x :: Real)
    xx = Float64(x)

    if xmin <= xx <= xmax

      return _maybe_nonneg(f_in(xx))

    elseif xx < xmin

      # -- left extrap
      xx >= xmin - Float64(Qpad_left) || error("left extrapolation detected beyond allowed window at x=$xx")

      y = if left_kind === :none
        error("No left extrapolation defined at x=$xx")
      elseif left_kind === :tangent
        yL + mL * (xx - xmin)
      elseif left_kind === :constant
        yL
      elseif left_kind === :secant
        yL + sL * (xx - xmin)
      elseif left_kind === :poly2_edge
        # -- use the interpolating polynomial if it's there
        f_in isa PolyInterp1D ? f_in(xx) :
          polyL === nothing ? yL + mL * (xx - xmin) : polyL(xx)
      elseif left_kind === :harmonic
        harmonic_model === nothing && error("No harmonic model available for left extrapolation")
        harmonic_model(xx)
      elseif left_kind === :zero
        0.0
      else
        error("Unsupported left extrapolation kind $left_kind")
      end

      return _maybe_nonneg(y)

    else

      # -- right extrap
      xx <= xmax + Float64(Qpad_right) || error("right extrapolation detected beyond allowed window at x=$xx")

      y = if right_kind === :none
        error("No right extrapolation defined at x=$xx")
      elseif right_kind === :tangent
        yR + mR * (xx - xmax)
      elseif right_kind === :constant
        yR
      elseif right_kind === :secant
        yR + sR * (xx - xmax)
      elseif right_kind === :poly2_edge
        f_in isa PolyInterp1D ? f_in(xx) :
          polyR === nothing ? yR + mR * (xx - xmax) : polyR(xx)
      elseif right_kind === :harmonic
        harmonic_model === nothing && error("No harmonic model available for right extrapolation")
        harmonic_model(xx)
      elseif right_kind === :zero
        0.0
      else
        error("Unsupported right extrapolation kind $right_kind")
      end

      return _maybe_nonneg(y)

    end

  end

  return fout

end

@inline _mode_V_override(
      mode_V_overrides :: AbstractDict{String, <:Any}
    , mode_name :: AbstractString
) = get(mode_V_overrides, String(mode_name), nothing)

@inline _res_interp_override(res_interp_overrides, channel_name, mode_name) =
  get(res_interp_overrides, (String(channel_name), String(mode_name)), nothing)

@inline _res_extrap_override(res_extrap_overrides, channel_name, mode_name) =
  get(res_extrap_overrides, (String(channel_name), String(mode_name)), nothing)

  ###################
  ###################
  ### P U B L I C ###
  ##VVVVVVVVVVVVVVV##
  #vvvvvvvvvvvvvvvvv#

  """
      build_dissociation_path_autostart(surface, space; q0_candidates = nothing, range_pad = 0.0, grad_tol = 1e-8)

  Constructs a dissociatin path by trying several starting points `q0` until a non-empty path is found
  """
  function build_dissociation_path_autostart(
      surface :: ChannelSurface
    , spec :: DissCalcSpec
    ; q0_candidates :: Union{Nothing, AbstractVector} = nothing
    , range_pad :: Float64 = 0.0
    , grad_tol :: Float64 = 1e-8
    )

    cands = q0_candidates === nothing ? _get_default_q0_candidates(surface) : q0_candidates

    for q0 in cands
      path = build_dissociation_path(
          surface
        , spec
        ; q0=q0
        , range_pad=range_pad
        , grad_tol=grad_tol
     )
      isempty(path.s) || return path
    end

    # -- nothing worked; return old default
  return build_dissociation_path(
      surface
    , spec
    ; q0=zeros(length(surface.modes))
    , range_pad=range_pad
    , grad_tol=grad_tol
  )

end

"""
    build_dissociation_paths_autostart(surfaces, spec; kwargs..)

Build several dissociation paths with `build_dissociation_path_autostart()`
"""
function build_dissociation_paths_autostart(
      surfaces :: AbstractVector{<:ChannelSurface}
    , spec :: DissCalcSpec
    ; kwargs...
  )
  paths = DissociationPathway[]
  for surface in surfaces
    push!(paths, build_dissociation_path_autostart(surface, spec; kwargs...))
  end
  return paths
end

"""
    prepare_dissociation_model(widetable, spec; targets_by_mode, res_interp_overrides, res_extrap_overrides, mode_V_overrides, kwargs...)

Widetable -> dissociation model pipeline
"""
function prepare_dissociation_model(
    widetable            :: Table
  , spec                 :: DissCalcSpec
  ; targets_by_mode      :: AbstractDict{String, <:Any} = Dict{String, Any}()
  , res_interp_overrides :: AbstractDict{Tuple{String, String}, <:Any} = Dict{Tuple{String, String}, Any}()
  , res_extrap_overrides :: AbstractDict{Tuple{String, String}, <:Any} = Dict{Tuple{String, String}, Any}()
  , mode_V_overrides     :: AbstractDict{String, <:Any} = Dict{String, Any}()
  , path_autostart       :: Bool = false
  , kwargs...
  )
  validate_calc_spec(widetable, spec)
  longtable = widetable_to_longtable(widetable, spec)
  curves    = parse_mode_curves(longtable; qcol = spec.qcol)
  fits      = fit_mode_curves(
      curves
    , spec
    ; targets_by_mode = targets_by_mode
    , res_interp_overrides = res_interp_overrides
    , res_extrap_overrides = res_extrap_overrides
    , mode_V_overrides = mode_V_overrides
  )
  surfaces  = build_channel_surfaces(fits, spec)
  paths     = if path_autostart
    build_dissociation_paths(surfaces, spec; kwargs...)
  else
    build_dissociation_paths_autostart(surfaces, spec; kwargs...)
  end
  pathfits  = fit_paths(
      paths
    ; kind = spec.interp.path_kind
  )
  return (
      longtable = longtable
    , curves = curves
    , fits = fits
    , surfaces = surfaces
    , paths = paths
    , pathfits = pathfits
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
  ds, smax = spec.path.ds, spec.path.smax
  ds   > 0    || error("ds ($ds) must be > 0")
  smax > 0    || error("smax ($smax) must be > 0")
  ds   < smax || error("smax ($smax) must be > ds ($ds)")

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
    fit_mode_curve(curve; interp, extrap, target, mode_freqs)

Builds a `ModeFit` from a `ModeCurve`, constructing callable fits for the following :
- `Eres(Q)`
- `Γ(Q)`
- `V(Q)`

Arguments
---------
- `curve` : Resonance curve extracted from the longtable
- `interp` : `InterpSpec` describing how to fit/interpolate `Eres`, `Γ`, and `V`
- `extrap` : `ExtrapSpec` describing how to fit/extrapolate `Eres`, `Γ`, and `V`
- `target` : Optional `TargetState` for this mode; required when `interp.V_kind` needs the target potential curve
- `mode_freqs` : Optional mode frequencies; needed when `interp.V_kind` requries them, e.g., `interp_Vkind === :harmonic`
"""
function fit_mode_curve(
      curve       :: ModeCurve
    ; interp      :: InterpSpec
    , extrap      :: ExtrapSpec
    , target      :: Union{TargetState, Nothing} = nothing
    , mode_freqs  :: AbstractDict{String, <:Real} = Dict{String, Float64}()
    , nedge_left  :: Int = 5
    , nedge_right :: Int = 5
  )

  isempty(curve.Q) &&
    error("Can't fit empty ModeCurve for channel '$(curve.channel_name)', 'mode = $(curve.mode_name)'")

  let allowed_Γ = (:constant, :tangent, :secant, :zero)
    extrap.Γ_left_kind  in allowed_Γ || error("Unsupported Γ_left_kind  = $(extrap.Γ_left_kind)")
    extrap.Γ_right_kind in allowed_Γ || error("Unsupported Γ_right_kind = $(extrap.Γ_right_kind)")
  end

  # -- min and max Q defined by the initial non-extrapolated curve
  Eres_qmin, Eres_qmax = extrema(curve.Q)

  # ----------
  # -- Eres(Q)
  # ----------

  # -- interpolation
  _Eres_fit = _make_interp(curve.Q, curve.Eres; kind=interp.Eres_kind)

  # -- extrapolation
  Eres_fit =
    (extrap.Eres_left_kind === :none && extrap.Eres_right_kind === :none) ?
    _Eres_fit :
    _wrap_with_extrapolation(
        _Eres_fit
      , curve.Q
      , curve.Eres
      ; left_kind   = extrap.Eres_left_kind
      , right_kind  = extrap.Eres_right_kind
      , Qpad_left   = extrap.EresΓ_Qpad_left
      , Qpad_right  = extrap.EresΓ_Qpad_right
      , nedge_left  = nedge_left
      , nedge_right = nedge_right
    )

  # -------
  # -- Γ(Q)
  # -------

  # -- interpolation
  _Γ_fit = _make_interp(curve.Q, curve.Γ; kind=interp.Γ_kind)

  # -- extrapolation
  Γ_fit = _wrap_with_extrapolation(
        _Γ_fit
      , curve.Q
      , curve.Γ
      ; left_kind   = extrap.Γ_left_kind
      , right_kind  = extrap.Γ_right_kind
      , Qpad_left   = extrap.EresΓ_Qpad_left
      , Qpad_right  = extrap.EresΓ_Qpad_right
      , nedge_left  = nedge_left
      , nedge_right = nedge_right
      , clamp_nonnegative = true # how would you interpret Γ<0 ?
    )


  # -------
  # -- V(Q)
  # -------
  if interp.V_kind === :harmonic

    haskey(mode_freqs, curve.mode_name) || error("Need mode frequency for mode '$(curve.mode_name)' when V_kind = :harmonic")
    ω = Float64(mode_freqs[curve.mode_name])

    Qeq = if target === nothing
      0.0
    else
      target.Q[argmin(target.V)]
    end

    # -- dimensionless coordinate ! V(Q) = ωQ²/2, NOT ω²Q²/2
    Vref0 = 0.5 * ω * (0.0 - Qeq)^2
    harmonic_model = HarmonicPotential(ω, Qeq, Vref0)

    V_fit = harmonic_model

  elseif interp.V_kind in (:linear, :poly2, :pchip, :akima)

    target === nothing && error("Need a TargetState for mode '$(curve.mode_name)' when V_kind = $(interp.V_kind)")

    # -- linear interpolation to approximate s.t. V(0) = 0.0
    V0 = linterp(target.Q, target.V, 0.0)

    # -- interpolation
    _V_fit = _make_interp(target.Q, target.V .- V0; kind=interp.V_kind)

    # -- optional harmonic model for extrapolation
    if haskey(mode_freqs, curve.mode_name)
      ω = Float64(mode_freqs[curve.mode_name])
      Qeq = target.Q[argmin(target.V)]
      Vref0 = 0.5 * ω * (0.0 - Qeq)^2
      harmonic_model = HarmonicPotential(ω, Qeq, Vref0)
    else
      harmonic_model = nothing
    end

    Vref = target.V .- V0

    # -- extrapolation
    V_fit =
      (extrap.V_left_kind === :none && extrap.V_right_kind === :none) ?
      _V_fit :
      _wrap_with_extrapolation(
          _V_fit
        , target.Q
        , Vref
        ; left_kind   = extrap.V_left_kind
        , right_kind  = extrap.V_right_kind
        , Qpad_left   = extrap.V_Qpad_left
        , Qpad_right  = extrap.V_Qpad_right
        , nedge_left  = nedge_left
        , nedge_right = nedge_right
        , harmonic_model = harmonic_model
      )

  else
    error("Unsupported V_kind '$(interp.V_kind)'")
  end

  # -- determine effective usable Q-domain post extrapolation
  Eres_qmin = extrap.Eres_left_kind  === :none ? Eres_qmin : Eres_qmin - extrap.EresΓ_Qpad_left
  Eres_qmax = extrap.Eres_right_kind === :none ? Eres_qmax : Eres_qmax + extrap.EresΓ_Qpad_right

  if interp.V_kind === :harmonic
    V_qmin, V_qmax = (-Inf, Inf)
  else
    V_qmin, V_qmax = extrema(target.Q)
    V_qmin = extrap.V_left_kind  === :none ? V_qmin : V_qmin - extrap.V_Qpad_left
    V_qmax = extrap.V_right_kind === :none ? V_qmax : V_qmax + extrap.V_Qpad_right
  end

  Qmin = max(Eres_qmin, V_qmin)
  Qmax = min(Eres_qmax, V_qmax)

  Qmin < Qmax || error("""
    No overlapping Q-grid found for
      channel: '$(curve.channel_name)'
      mode: '$(curve.mode_name)'
      Eres domain: [$Eres_qmin, $Eres_qmax]
      V domain: [$V_qmin, $V_qmax]
  """)


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
    fit_mode_curves(curves, spec; targets_by_mode, surface_interp_overrides, surface_extrap_overrides)

Fit a whole vector of `ModeCurve`s into a vector of `ModeFit`s.
`targets_by_mode` should map mode_name => TargetState (because V(Q) in this sense
is a property of the normal mode, not of the channel)
"""
function fit_mode_curves(
    curves :: AbstractVector{<:ModeCurve}
  , spec :: DissCalcSpec
  ; targets_by_mode :: AbstractDict{String, <:TargetState} = Dict{String, TargetState}()
  , res_interp_overrides :: AbstractDict{Tuple{String, String}, <:Any} = Dict{Tuple{String, String}, Any}()
  , res_extrap_overrides :: AbstractDict{Tuple{String, String}, <:Any} = Dict{Tuple{String, String}, Any}()
  , mode_V_overrides :: AbstractDict{String, <:Any} = Dict{String, Any}()
  )

  fits = ModeFit[]

  for curve in curves
    interp_over = _res_interp_override(res_interp_overrides, curve.channel_name, curve.mode_name)
    extrap_over = _res_extrap_override(res_extrap_overrides, curve.channel_name, curve.mode_name)
    V_over      = _mode_V_override(mode_V_overrides, curve.mode_name)
    interp = _merge_res_interp_spec(spec.interp, interp_over)
    extrap = _merge_res_extrap_spec(spec.extrap, extrap_over)
    if V_over !== nothing
      interp = InterpSpec(
          interp.Eres_kind
        , interp.Γ_kind
        , isnothing(V_over.V_kind) ? interp.V_kind : V_over.V_kind
        , interp.path_kind
      )
      extrap = ExtrapSpec(
          extrap.Eres_left_kind
        , extrap.Eres_right_kind
        , extrap.Γ_left_kind
        , extrap.Γ_right_kind
        , extrap.EresΓ_Qpad_left
        , extrap.EresΓ_Qpad_right
        , isnothing(V_over.V_left_kind) ? extrap.V_left_kind : V_over.V_left_kind
        , isnothing(V_over.V_right_kind) ? extrap.V_right_kind : V_over.V_right_kind
        , isnothing(V_over.V_Qpad_left) ? extrap.V_Qpad_left : V_over.V_Qpad_left
        , isnothing(V_over.V_Qpad_right) ? extrap.V_Qpad_right : V_over.V_Qpad_right
      )
    end
    push!(fits,
      fit_mode_curve(
          curve
        ; interp = interp
        , extrap = extrap
        , target = get(targets_by_mode, curve.mode_name, nothing)
        , mode_freqs = spec.mode_freqs
      )
    )
  end

  return fits

end

"""
   fits_for_channel(fits, channel)

Extracts the `ModeFit`s belonging to the channel `channel`, keyed by mode name.
Matches on channel name, mode name, and label
"""
function fits_for_channel(
    fits :: AbstractVector{<:ModeFit}
  , channel :: Channel
  )
  out = Vector{eltype(fits)}(undef, length(channel.modes))
  for (i, cm) in pairs(channel.modes)
    idx = findfirst(f ->
      f.curve.channel_name == channel.name &&
      f.curve.mode_name == cm.mode_name &&
      f.curve.label == cm.label
    , fits
    )
    idx === nothing && error("No ModeFit found for channel='$(channel.name)', mode='$(cm.mode_name)', label='$(cm.label)'")
    out[i] = fits[idx]
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
    fits :: AbstractVector{<:ModeFit}
  , channel :: Channel
  , spec :: DissCalcSpec
  )

  fitvec = fits_for_channel(fits, channel)
  modes   = _surface_modes(channel)

  isempty(modes) && error("Channel '$(channel.name)' has no modes !")

  # -- reference values at Q=0, averaged over mode fits
  E0arr = Float64[_fit_Eres0(fit) for fit in fitvec]
  Γ0arr = Float64[_fit_Γ0(fit)    for fit in fitvec]

  Eres0 = mean(E0arr)
  Γ0 = spec.combine_Γ === :additive ? mean(Γ0arr) :
       spec.combine_Γ === :max ? maximum(Γ0arr) :
       throw(ArgumentError("Unsupported combine_Γ = $(spec.combine_Γ)"))

  return ChannelSurface(
      channel.name
    , modes
    , fitvec
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
    fits :: AbstractVector{<:ModeFit}
  , spec :: DissCalcSpec
  )
  return [build_channel_surface(fits, ch, spec) for ch in spec.channels]
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
    , q :: AbstractVector{<:Real}
  )

  nq, nfits = length(q), length(surface.fits)
  nq == nfits || throw(ArgumentError("nq ($nq) != nfits ($nfits)"))
  val = surface.Eres0
  @inbounds for i in eachindex(surface.fits)
    fit = surface.fits[i]
    qi = Float64(q[i])
    val += fit.Eres_fit(qi) - _fit_Eres0(fit)
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
    , q :: AbstractVector{<:Real}
  )
  nq, nfits = length(q), length(surface.fits)
  nq == nfits || throw(ArgumentError("nq ($nq) != nfits ($nfits)"))
  if surface.combine_Γ === :additive
    val = surface.Γ0
    @inbounds for i in eachindex(surface.fits)
      fit = surface.fits[i]
      qi  = Float64(q[i])
      val += fit.Γ_fit(qi) - _fit_Γ0(fit)
    end
  elseif surface.combine_Γ === :max
    val = -Inf
    @inbounds for i in eachindex(surface.fits)
      fit = surface.fits[i]
      qi = Float64(q[i])
      Γi = fit.Γ_fit(qi)
      val = max(val, Γi)
    end
  else
    throw(ArgumentError("Unsupported combine_Γ = $(surface.combine_Γ)"))
  end
  return val
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
  nq, nfits = length(q), length(surface.fits)
  nq == nfits || throw(ArgumentError("nq ($nq) != nfits ($nfits)"))
  val = 0.0
  @inbounds for i in eachindex(surface.fits)
    fit = surface.fits[i]
    qi  = Float64(q[i])
    val += fit.V_fit(qi)
    # val += fit.V_fit(qi) - _fit_V0(fit)
  end
  return val
end

"""
    U(surface, q)

Evaluates the dissociative surface U(q) = V(q) + Eres(q)
"""
@inline function U(surface :: ChannelSurface, q :: AbstractVector)
  return V(surface, q) + Eres(surface, q)
end

#################
### GRADIENTS ###
#################

"""
    ∇E(surface, q)

Computes the gradient of Eres with respect to channel coordinates
"""
function ∇Eres(
    surface :: ChannelSurface
  , q :: AbstractVector
  )
  nq, nfits = length(q), length(surface.fits)
  nq == nfits || throw(ArgumentError("nq ($nq) != nfits ($nfits)"))
  g = Vector{Float64}(undef, nfits)
  @inbounds for i in eachindex(surface.fits)
    fit = surface.fits[i]
    qi  = Float64(q[i])
    g[i] = _deriv(fit.Eres_fit, qi)
  end
  return g
end

"""
    ∇V(surface, q)

Computes the gradient of V with respect to channel coordinate
"""
function ∇V(
    surface :: ChannelSurface
  , q :: AbstractVector
  )
  nq, nfits = length(q), length(surface.fits)
  nq == nfits || throw(ArgumentError("nq ($nq) != nfits ($nfits)"))
  g = Vector{Float64}(undef, nfits)
  @inbounds for i in eachindex(surface.fits)
    fit = surface.fits[i]
    qi  = Float64(q[i])
    g[i] = _deriv(fit.V_fit, qi)
  end
  return g
end

"""
    ∇U(surface, q)

Computes the gradient of U = V + Eres with respect to channel coordinate
"""
@inline function ∇U(surface :: ChannelSurface, q :: AbstractVector)
  return ∇V(surface, q) + ∇Eres(surface, q)
end

"""
   dUds(surface, q)

Calculates the derivative of U with respect to the steepest descent coordinate s:
   dq/ds = -∇U/|∇U| (vector quantity because q is a vector)
   dU/ds = ∇U · dq/ds = -|∇U|
"""
@inline function dUds(surface :: ChannelSurface, q :: AbstractVector)
  return -norm(∇U(surface, q))
end

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
  for i in eachindex(surface.fits)
    fit = surface.fits[i]
    (qv[i] < fit.qmin - pad || qv[i] > fit.qmax + pad) && return false
  end
  return true
end

"""
    _build_dissociation_path_partial(surfacde, spec; q0, direction=1, range_pad=0.0, grad_tol=1e-8)

Build the dissociation path for s≥0 (direction=1) or s≤0 (direction=-1)
"""
function _build_dissociation_path_partial(
    surface :: ChannelSurface
  , spec :: DissCalcSpec
  ; q0 :: AbstractVector = zeros(length(surface.modes))
  , direction :: Int = 1 # 1 -> s>0; -1 -> s<0
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

    q    .-= direction * spec.path.ds .* g ./ gnorm
    sval  += direction * spec.path.ds

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
    build_dissociation_path(surface, spec; q0, range_pad = 0.0, grad_tol=1e-8)

Builds the dissociation paths for all available s by merging the s≥0 and s≤0 paths.
"""
function build_dissociation_path(
    surface :: ChannelSurface
  , spec :: DissCalcSpec
  ; q0 :: AbstractVector = zeros(length(surface.modes))
  , range_pad :: Float64 = 0.0
  , grad_tol :: Float64 = 1e-8
  )
  pos = _build_dissociation_path_partial(
    surface, spec; q0=q0, direction= 1, range_pad=range_pad, grad_tol=grad_tol
   )
  neg = _build_dissociation_path_partial(
    surface, spec; q0=q0, direction=-1, range_pad=range_pad, grad_tol=grad_tol
   )
  return _merge_partial_dissociation_paths(neg, pos)
end

"""
    build_dissociation_paths(surface, spec; kwargs...)

Wrapper for multiple dissociation pathways.
"""
function build_dissociation_paths(
    surfaces :: AbstractVector{<:ChannelSurface}
  , spec :: DissCalcSpec
  ; kwargs...
  )
  pathways = DissociationPathway[]
  for surface in surfaces
    push!(pathways, build_dissociation_path(surface, spec; kwargs...))
  end
  return pathways
end

#######################
### PATH EVALUATION ###
### IN COORDINATE s ###
#######################

"""
    fit_path(path, kind)

Returns a DissociationPathwayFit; fits to all relevant quantities for a given pathway
"""
function fit_path(path :: DissociationPathway; kind :: Symbol=:linear)
  return DissociationPathwayFit(
    path
  , path.channel
  , path.s
  , path.q
  , _make_interp(path.s, path.Eres; kind = kind)
  , _make_interp(path.s, path.Γ;    kind = kind)
  , _make_interp(path.s, path.V;    kind = kind)
  , _make_interp(path.s, path.U;    kind = kind)
  , _make_interp(path.s, path.dUds; kind = kind)
  )
end

"""
    fit_paths(pathsl kind=:linear)

Returns an array of DissociationPathwayFit
"""
function fit_paths(
      paths :: AbstractVector{<:DissociationPathway}
    ; kind :: Symbol = :linear
  )
  res = DissociationPathwayFit[]
  for path in paths
    push!(res, fit_path(path; kind=kind))
  end
  return res
end


"Fitted Eres(s)"
Eres(pathfit :: DissociationPathwayFit, s :: Real) = pathfit.Eres_fit(_require_s_in_range(pathfit, s))
"Fitted Γ(s)"
Γ(pathfit    :: DissociationPathwayFit, s :: Real) = pathfit.Γ_fit(_require_s_in_range(pathfit, s))
"Fitted V(s)"
V(pathfit    :: DissociationPathwayFit, s :: Real) = pathfit.V_fit(_require_s_in_range(pathfit, s))
"Fitted U(s)"
U(pathfit    :: DissociationPathwayFit, s :: Real) = pathfit.U_fit(_require_s_in_range(pathfit, s))
"Fitted dUds(s)"
dUds(pathfit :: DissociationPathwayFit, s :: Real) = pathfit.dUds_fit(_require_s_in_range(pathfit, s))
