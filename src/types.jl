export EigenphData, EigenphRecord
export Resonance, IndexedResonance
export Table

"""
    EigenphData

single eigenphase dataset read from a single file
- E: energy grid (ne-array)
- δ: eigenphases (ne × nδ)
- names: descriptor for each column of δ
"""
struct EigenphData
  E :: Vector{Float64}
  δ :: Matrix{Float64}
  names :: Vector{String}
end

"""
    EigenphRecord

Generic container, attaches arbitrary metadata to EigenphData (irrep, mode..)
"""
struct EigenphRecord{M}
  data :: EigenphData
  path :: String
  meta :: M
end

"""
    Resonance

A resonance (candidate) extracted from EigenphData
- E: resonance position
- Γ: resonance width estimate (full width half max from |dδ/dE|)
- peak: peak heigh of |dδ/dE| for estimating strength
- col: eigenphase column index (data.names[col])
"""
struct Resonance
  E    :: Float64
  Γ    :: Union{Float64, Missing}
  peak :: Float64
  col  :: Int  # data.names[col]
end

"""
    Table

Minimal table container:
- header: Vector{String}
- rows: Matrix{Any}
- colmap: Dict{String, Int}
"""
struct Table
  header :: Vector{String}
  rows   :: Matrix{Any}
  colmap :: Dict{String, Int}
end

"""
    IndexedResonance

Just a Resonance with its index
"""
const IndexedResonance = Tuple{Int, Resonance}
