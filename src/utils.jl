using Printf

export select_cols, resolve_cols, resolve_resonances
export resonances_by_record
export write_table
export resonance_table_wide, resonance_table_long
export meta
export get_modenames, load_mode, load_modes
export width2lifetime

###############
### HELPERS ###
###############

"""
    meta(rec, key)

Helper forgetting metadata, also useful as an export
"""
meta(rec :: EigenphRecord, key :: Symbol; default = missing) = hasproperty(rec.meta, key) ? getproperty(rec.meta, key) : default

"""
    _yseries(data,  col; ytype, smoothwin)

returns the appropriate ytype series given Eigenphase data
- ytype :: Symbol
  - :phase (δ)
  - :deriv (dδ/dE)
  - :absderiv (|dδ/dE|)
- smoothwin :: Int -> moving average smoothing window
"""
function _yseries(
    data      :: EigenphData
  , col       :: Int
  ; ytype     :: Symbol = :phase
  , smoothwin :: Int  = 1
  )
  δ = @view data.δ[:,col]
  ytype === :phase && return δ
  dδ = _dydE(data.E, δ)
  dδ = smoothwin > 1 ? _moving_average(dδ, smoothwin) : dδ
  ytype === :deriv && return dδ
  ytype === :absderiv && return abs.(dδ)
  error("Symbol ytype must be :phase, :deriv, or :absderiv")
end

"Return non-missing (Q,y) sorted by Q"
function _clean_xy(Q :: AbstractVector, y :: AbstractVector)
  keep = .!ismissing.(y)
  Q2, y2 = Float64.((Q[keep], y[keep]))
  idx = sortperm(Q2)
  return Q2[idx], y2[idx]
end

"Converts units of energy E"
@inline function _convertE(E::Real, from_unit::Symbol, to_unit::Symbol)
  from_unit in ALLOWED_EUNITS || error("from_unit $from_unit is not in the list of allowed energy uints: $FROM_EUNITS")
  to_unit   in ALLOWED_EUNITS || error("to_unit $to_unit is not in the list of allowed energy uints: $FROM_EUNITS")
  from_unit === to_unit && return E
  from_unit === :eV && (E /= AU2EV) # FROM -> au
  to_unit   === :eV && (E *= AU2EV) # au -> TO
  return E
end
@inline _convertInvE(invE :: Real, from_unit :: Symbol, to_unit :: Symbol) =
  1/_convertE(invE, from_unit, to_unit)

##############
### PUBLIC ###
##############

# -- select columns by their names in various ways
# selector is a function
select_cols(data :: EigenphData, selector :: Function) =
  findall(selector, data.names)
# selector is a string
select_cols(data :: EigenphData, selector :: AbstractString) =
  findall(name -> occursin(selector, name), data.names)
# selector is an array of strings
select_cols(data :: EigenphData, selector :: AbstractVector{<:AbstractString}) =
  findall(name -> any(s-> occursin(s, name), selector), data.names)
# selector is a regex
select_cols(data :: EigenphData, selector :: Regex) =
  findall(name -> occursin(selector, name), data.names)

# convenience: selector = varargs strings, e.g. select_cols(data, "singlet A1", "triplet B2")
select_cols(data :: EigenphData, needles :: AbstractString...) =
    select_cols(data, collect(needles))

"""
    resolve_cols(data :: EigenphData; cols=:all, selector=nothing)

Resolve the columns from data.
- cols (:all or vector of indices)
- selector (string or selector function [name -> Bool])
"""
function resolve_cols(
    data :: EigenphData
  ; cols=:all
  , selector :: Union{Nothing, AbstractString, Function, Regex, AbstractVector{<:AbstractString}} = nothing
  )
  selector === nothing && return cols === :all ? collect(1:size(data.δ, 2)) : collect(cols)
  selector isa Union{AbstractString, Function, AbstractVector{<:AbstractString}, Regex} && return select_cols(data, selector)
  error("selector must be a Nothing, String, Regex, or a Function")
end

function resolve_resonances(data :: EigenphData, resonances, colidx :: AbstractVector{<:Integer}, res_detect :: NamedTuple)
  resonances === nothing && return nothing
  if resonances === :auto
    return find_resonances(data; cols=colidx, selector=nothing, res_detect...)
  elseif resonances isa AbstractVector{<:Resonance}
    return resonances
  else
    error("resonances must be nothing, :auto, or a vector of Resonance")
  end
end

"""
    resonances_by_record(idxres, nrecords)

Group IndexedResonance = (recid, Resonance) into per-record resonance lists.
Each list is sorted by resonance energy.
"""
function resonances_by_record(idxres :: AbstractVector{<:IndexedResonance}, nrecords :: Int)
  # -- array of arrays
  out = [Resonance[] for _ in 1:nrecords]
  for (recid, r) in idxres
    # -- skip any issues
    1 <= recid <= nrecords || continue
    push!(out[recid], r)
  end
  # -- sort by energy
  for o in out
    sort!(o, by = x -> x.E)
  end
  return out
end

# -- symlinks and directories under data_root will be taken to be normal modes
get_modenames(data_root :: AbstractString) = sort([d for d in readdir(data_root) if isdir(joinpath(data_root, d))])

"""
    load_mode(data_root, mode; sortQ, preprocess, Emin_unwrap, Eref_center)

Returns the records found in data_root/mode, where mode is a directory containing
a UKRmol+ run for a normal mode that lives inder the directory data_root.
"""
function load_mode(
        data_root   :: AbstractString
      , mode        :: AbstractString
      ; sortQ       :: Bool = true
      , preprocess  :: Bool = true
      , Emin_unwrap :: Real = 0.01
      , Eref_center :: Real = 0.01
      , period      :: Real = pi
    )

    root = joinpath(data_root, mode)
    isdir(root) || error("Mode directory DNE: $root")
    recs = DEAR.load_eigenph_records(root)

    sortQ && sort!(recs, by = r -> (ismissing(r.meta.Q) ? Inf : r.meta.Q))

    if preprocess
      DEAR.unwrap!(recs; Emin=Emin_unwrap, period=period)
      DEAR.center_records_at!(recs; Eref=Eref_center, period=period)
    end

    # -- get records, and attach mode into metadata
    recs = [DEAR.EigenphRecord(r.data, r.path, (; r.meta..., mode=mode)) for r in recs]

    return NormalMode(mode, root, recs)

end

"""
Wrapper for `load_mode` in case we just give  it a single directory that has the
dirname and modename, e.g., data/BEND
"""
function load_mode(path_or_root :: AbstractString; kwargs...)
  root = path_or_root
  isdir(root) || error("Mode directory  DNE: $root")
  mode = splitpath(normpath(root))[end]
  return load_mode(dirname(root), mode; kwargs...)
end

"Load all modes under data_root, and optionally pass a list of mode names"
function load_modes(data_root :: AbstractString; modes = :auto,  skip_missing :: Bool = false, kwargs...)
  names = modes === :auto ? get_modenames(data_root) : collect(modes)
  res = NormalMode[]
  for name in names
    root = joinpath(data_root, name)
    if !isdir(root)
      skip_missing && continue
      error("Mode directory DNE: $root")
    end
    push!(res, load_mode(data_root, name; kwargs...))
  end
  return res
end

"""
    width2lifetime(Γ; unit = :s)

Return τ=1/Γ (Inf if Γ=0, missing if Γ=missing).
`unit` should be one of
- :s (seconds)
- :fs (femtoseconds)
- :au (atomic units)
"""
function width2lifetime(Γ :: Union{Missing, Real}; unit :: Symbol = :s)
  ismissing(Γ) && return missing
  Γ == 0.0 && return Inf
  τ = 1 / Float64(Γ)
  return if unit === :au
    τ
  elseif unit === :s
    τ * AU2SEC
  elseif unit === :fs
    τ * AU2SEC*1e15
  else
    @error("Unexpected unit $unit; try :s, :fs, or :au")
  end
end
