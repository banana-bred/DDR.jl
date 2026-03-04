export EigenphData, EigenphRecord

"""
    EigenphData

single eigenphase dataset read from a single file
- E: energy grid (ne-array)
- δ: eigenphases (ne × nδ)
- name: descriptor for each column of δ
"""
struct EigenphData
  E :: Vector{Float64}
  δ :: Matrix{Float64}
  name :: Vector{String}
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
