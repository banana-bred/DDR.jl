using Statistics: median
export find_resonances, scan_resonances
###############
### HELPERS ###
###############

@inline function _range_indices(E :: AbstractVector; Emin=nothing, Emax=nothing)
  ne = length(E)
  i1 = Emin === nothing ?  1 : max(1,  searchsortedfirst(E, Emin))
  i2 = Emax === nothing ? ne : min(ne, searchsortedlast(E,  Emax))
  return i2 < i1 ? (1:0) : (i1:i2)
end

"""
    dydE(E, y)

Estimates the derivative of y(E) at E
"""
function dydE(E :: AbstractVector, y :: AbstractVector)
  ne = length(E)
  @assert length(y) == ne
  ne < 2 && return Float64[]

  dy = Vector{Float64}(undef, ne)
  # -- left
  dy[1] = (y[2] - y[1]) / (E[2] - E[1])
  # -- right
  dy[ne] = (y[ne] - y[ne-1]) / (E[ne] - E[ne-1])
  # -- the rest
  for i in 2:ne-1
    dy[i]= (y[i+1] - y[i-1]) / (E[i+1] - E[i-1])
  end
  return dy
end

"""
    _moving_average(y, w)

Moving average smoothing of y over an odd window of w values
"""
function _moving_average(y :: AbstractVector, w :: Int)
  w <= 1 && return copy(y)
  isodd(w) || error("Smoothing window w must be odd")
  n = length(y)
  out = Vector{Float64}(undef, n)
  hw = div(w-1, 2)
  for i in 1:n
    lo = max(1, i - hw)
    hi = min(n, i + hw)
    out[i] = sum(@view y[lo:hi]) / (hi - lo + 1)
  end
  return out
end

"""
    _get_local_maxima_idx(y)

Returns the indices of local maxima in y
"""
function _get_local_maxima_idx(y :: AbstractVector)
  n = length(y)
  idx = Int[]
  n < 3 && return idx
  for i in 2:n-1
    y[i] > y[i-1] && y[i] >= y[i+1] && push!(idx, i)
  end
  return idx
end

"""
    _prominence(y, ip, win)

Estimates the prominence (y[ip] - max(data_in_window)) of data in y at index ip with respect to the window ±win
"""
function _prominence(y :: AbstractVector, ip :: Int, win :: Int)
  n = length(y)
  lo = max(1, ip - win)
  hi = min(n, ip + win)
  lmin = minimum(@view y[lo:ip])
  rmin = minimum(@view y[ip:hi])
  return y[ip] - max(lmin, rmin)
end

"""
    _enforce_min_peak_spacing(E, peaks, height; min_distance)

Tries to enforce a minimum spacing between resonances
"""
function _enforce_min_peak_spacing(
      E            :: AbstractVector
    , peaks        :: Vector{Int}
    , height       :: AbstractVector
    ; min_distance :: Real = 0.0
  )
  min_distance <= 0 && return peaks
  # -- biggest resonances first
  ord = sort(peaks, by = i -> height[i], rev=true)
  kept = Int[]
  for i in ord
    keep = true
    # -- keep the first, but enforce spacing for the rest
    for j in kept
      if abs(E[i] - E[j]) < min_distance
        keep = false
        break
      end
    end
    keep && push!(kept, i)
  end
  sort!(kept)
  return kept
end

"""
    _fwhm(E, y, ip)

Returns the width of a resonance at y[ip]
"""
function _fwhm(E :: AbstractVector, y :: AbstractVector, ip :: Int)
  n = length(y)
  y[ip] <= 0 && return missing
  half = y[ip]/2

  # -- left half point
  ileft = nothing
  for i in ip-1:-1:1
    y[i] > half && continue
    ileft = i
    break
  end

  # -- right half point
  iright = nothing
  for i in ip+1:n
    y[i] > half && continue
    iright = i
    break
  end

  (ileft === nothing || iright === nothing) && return missing

  # -- linear interpolation left
  x1, x2 = E[ileft], E[ileft+1]
  y1, y2 = y[ileft], y[ileft+1]
  y1 == y2 && return missing
  t = (half-y1)/(y2-y1)
  Eleft = x1 + t*(x2-x1)

  # -- linear interpolation right
  x1, x2 = E[iright-1], E[iright]
  y1, y2 = y[iright-1], y[iright]
  y1 == y2 && return missing
  t = (half-y1)/(y2-y1)
  Eright = x1 + t*(x2-x1)

  return Eright - Eleft

end

##############
### PUBLIC ###
##############

"""
    function find_resonances(data :: EigenphData
      ; cols=:all
      , selector=nothing
      , Emin=nothing
      , Emax=nothing
      , smoothwin :: Int = 1
      , min_height :: Real = 0.0
      , min_prominence :: Real = 0.0
      , prominence_window :: Int = 25
      , min_distance :: Union{Real, Nothing} = 0.0
      )

Detect resonances using peaks in |dδ/dE|

Notes:
- `smoothwin` is an odd window in E for a moving average procedure on dδ/dE. Disabled when =1
"""
function find_resonances(data :: EigenphData
  ; cols=:all
  , selector=nothing
  , Emin=nothing
  , Emax=nothing
  , smoothwin :: Int =1
  , min_height :: Real = 0.0
  , min_prominence :: Real = 0.0
  , prominence_window :: Int = 25
  , min_distance :: Union{Real, Nothing} = nothing
  )

  res = Resonance[]

  colidx = resolve_cols(data; cols=cols, selector=selector)
  isempty(colidx) && return res

  E = data.E
  ie_range = _range_indices(E; Emin=Emin, Emax=Emax)

  isempty(ie_range) && return res

  for col in colidx
    δcol = @view data.δ[:, col]
    dδdE = dydE(E, δcol)
    dδdE = smoothwin > 1 ? _moving_average(dδdE,  smoothwin) : dδdE

    # -- filter based on energy subinterval
    Esub = @view E[ie_range]
    yabs = abs.(@view dδdE[ie_range])

    peaks = _get_local_maxima_idx(yabs)

    # -- filter peaks based on threshold + prominence
    peaks = [ip for ip in peaks if yabs[ip] >= min_height
      && _prominence(yabs, ip, prominence_window) >= min_prominence]

    # -- minimum energy spacing
    mindist = if min_distance === nothing
      ΔE = median(diff(Esub))
      max(8, 2*smoothwin)*ΔE
    else
      min_distance
    end
    peaks = _enforce_min_peak_spacing(Esub, peaks, yabs; min_distance=mindist)

    for ip in peaks
      # -- map peak index ip to global energy index
      ig = first(ie_range) - 1 + ip
      Eres = E[ig]
      peak = yabs[ip]
      Γ = _fwhm(Esub, yabs, ip)
      push!(res,Resonance(Eres, Γ, peak, col))
    end

  end

  sort!(res, by = r -> r.E)
  return res

end

"""
    function scan_resonances(records :: AbstractVector{<:EigenphRecord}; kwargs...)

Returns (recid, resonance) pairs where `recid` indexes `records`, so that we know
where this resonance comes from without duplicating too much data.
"""
function scan_resonances(records :: AbstractVector{<:EigenphRecord}; kwargs...)
  hits = Tuple{Int, Resonance}[]
  for (recid, r) in pairs(records)
    for res in find_resonances(r.data; kwargs...)
      push!(hits, (recid,res))
    end
  end
  return hits
end
