export select_cols

# -- select columns by their names
select_cols(data :: EigenphData, selector :: Function) =
  findall(selector, data.name)
select_cols(data :: EigenphData, needle :: AbstractString) =
  findall(name -> occursin(needle, name), data.name)
