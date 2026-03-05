export select_cols, resolve_cols

# -- select columns by their names
select_cols(data :: EigenphData, selector :: Function) =
  findall(selector, data.names)
select_cols(data :: EigenphData, needle :: AbstractString) =
  findall(name -> occursin(needle, name), data.names)

"""
    resolve_cols(data :: EigenphData; cols=:all, selector=nothing)

Resolve the columns from data.
- cols (:all or vector of indices)
- selector (string or selector function [name -> Bool])
"""
function resolve_cols(data :: EigenphData; cols=:all, selector=nothing)
  selector === nothing && return cols === :all ? collect(1:size(data.δ, 2)) : collect(cols)
  selector isa Union{AbstractString, Function} && return select_cols(data, selector)
  error("selector must be a String of a Function")
end
