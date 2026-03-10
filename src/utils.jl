using Printf

export select_cols, resolve_cols, resolve_resonances
export resonances_by_record
export write_table
export resonance_table_wide, resonance_table_long
export meta

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
  findall(name -> any(s-> occursin(s, name), selector))
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
  , selector :: Union{Nothing, AbstractString, Function, Regex} = nothing
  )
  selector === nothing && return cols === :all ? collect(1:size(data.δ, 2)) : collect(cols)
  selector isa Union{AbstractString, Function, Regex} && return select_cols(data, selector)
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
