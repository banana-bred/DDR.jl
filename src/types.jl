export EigenphData, EigenphRecord
export Resonance, IndexedResonance
export Table
export NormalMode
export Channel, ChannelMode
export TargetEnergies, TargetState
export DissociationPathway, DissociationPathwayFit
export HarmonicPotential

"""
    EigenphData

single eigenphase dataset read from a single file
- E: energy grid (ne-array)
- δ: eigenphases (ne × nδ)
- names: descriptor for each column of δ
"""
struct EigenphData
  E     :: Vector{Float64}
  δ     :: Matrix{Float64}
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
  name    :: String
  root    :: String
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

"""
    ChannelMode

The name of a mode to which a channel is attached, with a label in case there are several,
e.g., BEND_1, BEND_2.
"""
struct ChannelMode
  mode_name :: String
  label     :: String
end

"""
    Channel

One of the channels in the calculation, which can be attached to several channels
"""
struct Channel
  name  :: String
  modes :: Vector{ChannelMode}   # -- enabled modes in coordinate order
end

"""
    InterpSpec

Interpolation specification. The kinds of interpolation that are performed
for the different compontents of the approach.
"""
struct InterpSpec
  Eres_kind :: Symbol # -- :linear, :poly2, :spline
  Γ_kind    :: Symbol # -- :linear, :poly2, :spline
  V_kind    :: Symbol # -- :harmonic, :linear, :poly2, :spline
  path_kind :: Symbol # -- :linear, :poly2..
end

"""
    PathSpec

Dissociation pathway specification, where `s` is the dissociation coordinate.
- `ds`: step size
- `smax`: larges value for `abs(s)` considered
"""
struct PathSpec
  ds   :: Float64 # -- path step size
  smax :: Float64 # -- maxi path length
end

"""
    DissCalcSpec

Dissociation calculation specification
"""
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

"""
    ModeCurve

"""
struct ModeCurve
  channel_name :: String
  mode_name    :: String
  label        :: String
  Q            :: Vector{Float64}
  Eres         :: Vector{Float64}
  Γ            :: Vector{Float64}
end

"""
    ModeFit
"""
struct ModeFit{FE, FΓ, FV}
  curve    :: ModeCurve
  target   :: Union{TargetState, Nothing}
  Eres_fit :: FE
  Γ_fit    :: FΓ
  V_fit    :: FV
  qmin     :: Float64
  qmax     :: Float64
end

"""
    ChannelSurface

The surface of a particular channel; combined resonances across modes.

For example, a surface can be made of
a 3B₂ resonance in the BEND mode and a 3B₂ resonance in the SYMSTRETCH mode.
Another surface may have the first 3B₂ resonance in the BEND mode and another 3B₂ resonance
in the SYMSTRETCH mode. Yet another surface may have only the 1A₁ resonance in the BEND mode.
"""
struct ChannelSurface{MF<:ModeFit}
  channel_name :: String
  modes        :: Vector{String}        # -- defines the coordinate order
  fits         :: Vector{MF} # same order as modes
  mode_freqs   :: Dict{String, Float64}
  combine_Eres :: Symbol
  combine_Γ    :: Symbol
  Eres0        :: Float64
  Γ0           :: Float64
end

"""
    DissociationPathway

The 1D dissociation pathway corresponding to the steepest descent along a ChannelSurface
"""
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

"""
    DissociationPathwayFit

The fits to a DissociationPathway
"""
struct DissociationPathwayFit{FE, FΓ, FV, FU, FdU}
  path     :: DissociationPathway # original path
  channel  :: String
  s        :: Vector{Float64}
  q        :: Matrix{Float64}
  Eres_fit :: FE
  Γ_fit    :: FΓ
  V_fit    :: FV
  U_fit    :: FU
  dUds_fit :: FdU
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

"""
    LinearInterp1D

A 1D linear interpolation
"""
struct LinearInterp1D
  x :: Vector{Float64}
  y :: Vector{Float64}
end

"""
    PolyInterp1D

A 1D polynomila interpolation.
Coefficients are stored in the `coeffs` array, starting with power 0
"""
struct PolyInterp1D
  coeffs :: Vector{Float64}
  xmin   :: Float64
  xmax   :: Float64
  clampx :: Bool
end

struct HarmonicPotential
  ω     :: Float64
  Qeq   :: Float64
  Vref0 :: Float64
end
@inline function (h :: HarmonicPotential)(Q :: Float64)
  return 0.5*h.ω * (Q - h.Qeq)^2 - h.Vref0
end
@inline (h :: HarmonicPotential)(Q :: Real) = h(Float64(Q))
