using LinearAlgebra: \
using DataInterpolations: PCHIPInterpolation, AkimaInterpolation

###############
### HELPERS ###
###############

"""
    _sorted_xy(x, y; rev=false)

Return (x,y), sorted by `x` ascending (rev=false) or descending (rev=true)
"""
function _sorted_xy(x :: AbstractVector, y :: AbstractVector; rev :: Bool = false)
  length(x) == length(y) || error("x/y length mismatch")
  idx = sortperm(x; rev=rev)
  return Float64.(x[idx]), Float64.(y[idx])
end

"""
    _polyfit(x, y, deg)

Least-squares polynomial fit with coefficients in ascending powers:
    c₀ + c₁ x + c₂ x² + ..
"""
function _polyfit(x :: Vector{Float64}, y :: Vector{Float64}, deg :: Int)

  n = length(x)
  n >=2 || error("At least 2 points needed to fit a polynomial")
  d = min(deg, n-1)

  A = Matrix{Float64}(undef, n, d+1)
  A[:, 1] .= 1.0
  for j in 2:d+1
    @inbounds A[:, j] .= A[:, j-1] .* x
  end

  return A \ y

end

"""
    _make_interp(x, y; kind=:linear, clampx=false)

Constructs a simple 1D interpolant. Allowed `kind`s:
- :linear
- :poly2 (quadratic)
- :pchip
- :akima
- :spline (currently undefined/unsupported)
"""
function _make_interp(
    x :: AbstractVector
  , y :: AbstractVector
  ; kind :: Symbol =:linear
  , clampx :: Bool = false
  )

  xs, ys = _sorted_xy(x, y)

  if kind === :linear

    return LinearInterp1D(xs, ys)

  elseif kind === :poly2

    c = _polyfit(xs, ys, 2)
    return PolyInterp1D(c, xs[1], xs[end], clampx)

  elseif kind === :pchip

    itp = PCHIPInterpolation(ys, xs)
    return PCHIPInterp1D(itp, xs[1], xs[end], clampx)

  elseif kind === :akima

    itp = AkimaInterpolation(ys, xs)
    return AkimaInterp1D(itp, xs[1], xs[end], clampx)

  elseif kind === :spline

    error("Spline not programmed yet !")

  else

    error("Unsupported interpolation kind $kind")

  end

end
