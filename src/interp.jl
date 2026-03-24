using LinearAlgebra: \

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
    LinearInterp1D(xq)

Callable version of the LinearInterp1D struct
"""
(f::LinearInterp1D)(xq::Real) = linterp(f.x, f.y, Float64(xq))

"""
    PolyInterp1D(xq)

Callable version of the PolyInterp1D struct
"""
function (f :: PolyInterp1D)(xq :: Real)
  x = Float64(xq)
  f.clampx && (x = clamp(x, f.xmin, f.xmax))
  y = f.coeffs[end]
  for i in (length(f.coeffs) - 1):-1:1
    @inbounds y = muladd(y, x, f.coeffs[i])
  end
  return y
end

"""
    _make_interp(x, y; kind=:linear, clampx=true)

Constructs a simple 1D interpolant. Allowed `kind`s:
- :linear
- :poly2 (quadratic)
- :spline
"""
function _make_interp(
    x :: AbstractVector
  , y :: AbstractVector
  ; kind :: Symbol =:linear
  , clampx :: Bool =true
  )

  xs, ys = _sorted_xy(x, y)

  if kind === :linear
    return LinearInterp1D(xs, ys)
  elseif kind === :poly2
    c = _polyfit(xs, ys, 2)
    return PolyInterp1D(c, xs[1], xs[end], clampx)
  elseif kind === :spline
    error("Spline not programmed yet !")
  else
    error("Unsupported interpolation kind $kind")
  end
end
