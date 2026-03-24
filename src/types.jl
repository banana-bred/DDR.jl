export EigenphData, EigenphRecord
export Resonance, IndexedResonance
export Table
export NormalMode
export TargetEnergies, TargetState

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

struct NormalMode{M}
  name :: String
  root :: String
  records :: Vector{EigenphRecord{M}}
end

"""
    TargetEnergies

Contains all target state energies read from the file `target.energies` for each Q
"""
struct TargetEnergies
  Q     :: Vector{Float64}
  E     :: Matrix{Float64}
  names :: Vector{String}
  path  :: String
end

"""
    TargetState

A singlet target state potential energy curve V(Q)
"""
struct TargetState
  Q    :: Vector{Float64}
  V    :: Vector{Float64}
  name :: String
  col  :: Int
end

struct ChannelMode
  mode_name :: String
  label     :: String
end

"One of the channels to consider"
struct Channel
  name  :: String
  modes :: Vector{ChannelMode}   # -- enabled modes in coordinate order
end

"Interpolation specification"
struct InterpSpec
  Eres_kind :: Symbol # -- :linear, :poly2, :spline
  Γ_kind    :: Symbol # -- :linear, :poly2, :spline
  V_kind    :: Symbol # -- :harmonic, :linear, :poly2, :spline
end

"Dissociation path specification"
struct PathSpec
  ds   :: Float64 # -- path step size
  smax :: Float64 # -- maxi path length
end

"Dissociation calculation specification"
struct DissCalcSpec
  channels     :: Vector{Channel}
  mode_freqs   :: Dict{String, Float64}
  qcol         :: String
  energy_grid  :: Vector{Float64}
  combine_Eres :: Symbol # -- :additive, ..
  combine_Γ    :: Symbol # -- :additive, :max..
  interp       :: InterpSpec
  path         :: PathSpec
end

struct ModeCurve
  channel_name :: String
  mode_name    :: String
  label        :: String
  Q            :: Vector{Float64}
  Eres         :: Vector{Float64}
  Γ            :: Vector{Float64}
end

struct ModeFit
  curve  :: ModeCurve
  target :: Union{TargetState, Nothing}
  Eres_fit
  Γ_fit
  V_fit
  qmin   :: Float64
  qmax   :: Float64
end

struct ChannelSurface
  channel_name :: String
  modes        :: Vector{String}        # -- defines the coordinate order
  fits         :: Dict{String, ModeFit} # -- keyed by mode name
  mode_freqs   :: Dict{String, Float64}
  combine_Eres :: Symbol
  combine_Γ    :: Symbol
  Eres0        :: Float64
  Γ0           :: Float64
end

"Sampled dissociation pathway"
struct DissociationPathway
  channel :: String
  modes   :: Vector{String}  # -- matches columns of q
  s       :: Vector{Float64}
  q       :: Matrix{Float64} # -- (length(s), length(modes))
  Eres    :: Vector{Float64}
  Γ       :: Vector{Float64}
  V       :: Vector{Float64}
  U       :: Vector{Float64}
  dUds    :: Vector{Float64}
end

"Fitted dissociation pathway"
struct DissociationPathwayFit
  channel :: String
  s       :: Vector{Float64}
  q       :: Matrix{Float64}
  Eres
  Γ
  V
  U
  dUds
end

"Dissociation result from a subset of considered channels"
struct PartialDissociationResult
  channel  :: String
  energies :: Vector{Float64}
  σ        :: Vector{Float64}
  path     :: DissociationPathway
end

"Dissociation result from all considered channels, and their branching ratios"
struct TotalDissociationResult
  partials          :: Vector{PartialDissociationResult}
  energies          :: Vector{Float64}
  σ                 :: Vector{Float64}
  channel_fractions :: Dict{String, Vector{Float64}}
end

"A 1D linear interpolation"
struct LinearInterp1D
  x :: Vector{Float64}
  y :: Vector{Float64}
end

"A 1D polynomial interpolation"
struct PolyInterp1D
  coeffs :: Vector{Float64}
  xmin :: Float64
  xmax :: Float64
  clampx :: Bool
end
