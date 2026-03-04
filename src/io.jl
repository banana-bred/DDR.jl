using DelimitedFiles: readdlm

# -- procedures
export read_eigenph_file
export read_geom_file
export discover_eigenph_files
export load_eigenph_records
export geom_index
export select_cols

"""
    read_header_tokens(path; comment_char = '#', header_selector = tokens -> true)

Returns the header tokens from a file, like

    #Energy   singlet A1     triplet A1     singlet A2 .....
    #geom     Q

Optinally, uses `header_selector` to filter extra comment lines that may appear
to determine which line to use as the header in determining tokens.
"""
function _read_header_tokens(
      path :: AbstractString
    ; comment_char :: Char = '#'
    , header_selector = tokens -> true
  )
  return open(path, "r") do io
    for line in eachline(io)
      # -- skip emtpy lines
      s = strip(line)
      isempty(s) && continue
      # -- break when we find first non-comment line (allow for spaces before comment_char)
      t = lstrip(s)
      startswith(t, string(comment_char)) || break
      # -- strip comment char
      s2 = strip(chopprefix(t, string(comment_char)))
      isempty(s2) && continue
      tokens = split(s2)
      isempty(tokens) && continue
      header_selector(tokens) && return tokens
    end
  end
  return nothing
end

"""
    names_from_tokens(tokens :: Vector{String}, ncols  :: Int)

["singlet", "A1", "triplet", "A1"..] -> ["singlet A1", "triplet A1"..]
["Q"] -> ["Q"]
"""
function names_from_tokens(tokens :: AbstractVector{<:AbstractString}, ncols :: Int)
  ntokens =  length(tokens)
  if ntokens == ncols
    return tokens
  end
  if ntokens > 0 && ntokens %ncols == 0
    # -- integer division ntokens / ncols
    k = div(ntokens, ncols)
    return [join(tokens[1+(j-1)*k : j*k], " ") for j in 1:ncols]
  end
  return ["col$(j)" for j in 1:ncols]
end

"""
    read_eigenph_file(path :: AbstractString
      ; comments :: Bool = true
      , comment_char :: Char = '#'
      , header_selector = t -> lowercase(t[1]) == "energy"
    )

Returns the energies E and eigenphases δ from the tabular data at path
"""
function read_eigenph_file(
    path :: AbstractString
  ; comments :: Bool = true
  , comment_char :: Char= '#'
  , header_selector = t -> lowercase(t[1]) == "energy"
  )

  raw = readdlm(path, comments = comments, comment_char = comment_char)

  ndims(raw) == 2 || error("Expected a 2D table in $path, got $(ndims(raw))D")
  size(raw, 2) >= 2 || error("Expected ≥ 2 columns in $path, got $(size(raw, 2))")

  E = Float64.(raw[:, 1])
  δ = Float64.(raw[:, 2:end])
  nδ = size(δ,  2)

  tokens = _read_header_tokens(path
    ; comment_char = comment_char
    , header_selector = header_selector
  )

  if tokens === nothing
    names = ["phase$(iδ)" for iδ in 1:nδ]
  else
    # -- drop the first one
    names = names_from_tokens(tokens[2:end], nδ)
  end

  return EigenphData(E, δ, names)

end

"""
    read_geom_file(
        path :: AbstractString
      ; comments :: Bool = true
      , comment_char :: Char= '#'
      , header_selector = t -> occursin("geom", lowercase(t[1]))
    )

Returns the tuple (names,  raw):
- names from tokens of the geometry file using the header_selector function
- raw data obtained from file read while ignoring geometries
"""
function  read_geom_file(
    path :: AbstractString
  ; comments :: Bool = true
  , comment_char :: Char= '#'
  , header_selector = t -> occursin("geom", lowercase(t[1]))
  )
  raw = readdlm(path; comments=comments, comment_char=comment_char)
  tokens = _read_header_tokens(path; comment_char=comment_char, header_selector = header_selector)
  names = tokens === nothing ? ["cols$(iδ)" for iδ in 1:size(raw, 2)] : names_from_tokens(tokens, size(raw, 2))
  return names, raw
end

"""
    geom_value_map(raw; geom_col= 1, val_col= 2) -> Dict{Int, Float64}
"""
function geom_value_map(raw; geom_col :: Int = 1, val_col :: Int = 2)
  g = Int.(raw[:, geom_col])
  v = Float64.(raw[:, val_col])
  return Dict(g[i] => v[i] for i in eachindex(g))
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
      push!(files, joinpath(root, f))
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
