using DelimitedFiles: readdlm

# -- types
export EigenphData, EigenphRecord

# -- procedures
export read_eigenph_file
export discover_eigenph_files
export load_eigenph_records
export geom_index

"""
    EigenphData

single eigenphase dataset read from a single file
- E: energy grid (ne-array)
- δ: eigenphases (ne × nδ)
"""
struct EigenphData
  E :: Vector{Float64}
  δ :: Matrix{Float64}
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
    read_eigenph_file(path :: AbstractString; comments :: Bool = true, comment_char :: AbstractChar = '#')

Returns the energies E and eigenphases δ from the tabular data at path
"""
function read_eigenph_file(path :: AbstractString; comments :: Bool = true, comment_char :: AbstractChar= '#')

  raw = readdlm(path, comments = comments, comment_char = comment_char)

  ndims(raw) == 2 || error("Expected a 2D table in $path, got $(ndims(raw))D")
  size(raw, 2) >= 2 || error("Expected ≥ 2 columns in $path, got $(size(raw, 2))")

  E = Float64.(raw[:, 1])
  δ = Float64.(raw[:, 2:end])

  return EigenphData(E, δ)

end

"""
    geom_index(path :: AbstractString) -> Int || nothing

Extract the geometry number from a filename of the form ..geomNUM, where NUM is an integer.
If there's as match, return NUM. Otherwise, return nothing.
"""
function geom_index(path :: AbstractString)
  fname = splitpath(path)[end]
  m = match(r"geom(\d+)\s*$", fname)
  return m === nothing ? nothing : parse(Int, m.captures[1])
end

"""
    function discover_eigenph_files(root :: AbstractString ; recursive :: Bool = true, fname_regex :: Regex = r"eigenph\\.all\\.geom\\d+\$")

Find eigenphase files under the directory ROOT. Returns a vector of full file paths.
"""
function discover_eigenph_files(
    root        :: AbstractString
  ; recursive   :: Bool = true
  , fname_regex :: Regex = r"eigenph\.all\.geom\d+$"
  )


  files = String[]

  if recursive

    for (dir, _, fnames) in walkdir(root)
      for f in fnames
        occursin(fname_regex, f) || continue
        push!(files, joinpath(dir,  f))
      end
    end

  else

    for f in readdir(root)
      occursin(fname_regex, f) || continue
      push!(files, joinpath(dir, f))
    end

  end

  # -- sort by geometry index
  sort!(files, by = p -> something(geom_index(p), typemax(Int)))
  return files

end

"""
    load_eigenph_records(
        root :: AbstractString
      ; path_meta = (p -> (;))
      , recursive :: Bool = true
      , fnameRegex :: Regex = r"eigenph\\.all\\.geom\\d+\$"
      , comments :: Bool = true
    )

Discovers eigenphase files under `root`, reads them, and attaches appropriate metadata, using
the user-supplied function `path_meta(path)`.

Example for pathMeta: (mode="BEND", spinmult="singlet", irrep="A1")
"""
function load_eigenph_records(
      root :: AbstractString
    ; path_meta = (p -> (;))
    , recursive :: Bool = true
    , fname_regex :: Regex = r"eigenph\.all\.geom\d+$"
    , comments :: Bool = true
  )


  paths = discover_eigenph_files(root; recursive = recursive, fname_regex=fname_regex)
  records = Vector{EigenphRecord}(undef, length(paths))
  for (i, p) in pairs(paths)
    d = read_eigenph_file(p; comments=comments)
    meta = path_meta(p)
    records[i] = EigenphRecord(d, String(p), meta)
  end

  return records

end
