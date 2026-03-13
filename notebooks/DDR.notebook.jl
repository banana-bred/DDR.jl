### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9527d1fc-1f12-11f1-bbc3-ad70a776bfa7
begin
    import Pkg
    Pkg.activate(@__DIR__)                      # -- notebooks/Project.toml
    Pkg.develop(path=joinpath(@__DIR__, ".."))  # -- add local DDR package
    Pkg.instantiate()
    using DDR
	using Plots
	using PlutoUI
end

# ╔═╡ c47fdc16-11f3-4ac0-af04-f6468d48e989
md"# Direct Dissociative Recombination"

# ╔═╡ f6999a5f-9003-4270-ac1c-29a2350711ff
begin
	fname_regex = r"^eigenph\.all\.geom\d+$"

	function discover_modes(data_root::AbstractString;
	        require_eigenph::Bool=true,
	        recursive::Bool=false)
	
	    isdir(data_root) || return String[]
	
	    modes = String[]
	    for name in readdir(data_root)
	        startswith(name, ".") && continue
	        path = joinpath(data_root, name)
	        isdir(path) || continue
	
	        ok = true
	        if require_eigenph
	            ok = false
	
	            # fast, non-recursive check for your typical layout
	            eigenph_dir = joinpath(path, "collected_scattering_data", "eigenph")
	            if isdir(eigenph_dir)
	                ok = any(f -> occursin(fname_regex, f), readdir(eigenph_dir))
	            end
	
	            # optional deeper scan if layout varies
	            if !ok && recursive
	                for (dir, _, files) in walkdir(path; follow_symlinks=true)
	                    if any(f -> occursin(fname_regex, f), files)
	                        ok = true
	                        break
	                    end
	                end
	            end
	        end
	
	        ok && push!(modes, name)
	    end
	
	    sort!(modes)
	    return modes
	end
end

# ╔═╡ 94c35eb8-7e95-4f76-9829-94da62bc8878
md"""
### Data root dir: $(@bind data_root TextField(default="data"))
"""

# ╔═╡ cc01839f-ed8a-46e4-b112-8d91404ea61a

begin
mode_names = discover_modes(data_root; require_eigenph=true, recursive=false)
mode_options = isempty(mode_names) ? ["<no modes found>"] : mode_names

md"### Pick a normal mode: $(@bind mode Select(mode_options))"
end

# ╔═╡ 91911778-83ba-41c1-86a3-18ab1cb7c295
@bind cfg PlutoUI.combine() do Child
md""" ## DDR controls

**Selector** (regex string): $(Child(TextField(default="(singlet A1)")))

**ytype**: $(Child(Select(["phase","deriv","absderiv"])))

**Auto Emin**: $(Child(CheckBox(default=false)))
**Emin**: $(Child(Slider(0.0:0.01:2.0, default=0.05, show_value=true)))	

**Auto Emax**: $(Child(CheckBox(default=false)))
**Emax**: $(Child(Slider(0.0:0.01:10.0, default=1.5, show_value=true)))	

**min_height**: $(Child(Slider(0:1:200, default=40, show_value=true)))

**min_prominence**: $(Child(Slider(0:0.1:100, default=0.5, show_value=true)))

**smoothwin**: $(Child(Slider(1:2:21, default=1, show_value=true)))
"""
end

# ╔═╡ a0854639-7af8-42e9-95f4-9561b7a6e420
begin
	selector_re = Regex(cfg[1])                 
	ytype       = Symbol(cfg[2])                
	Emin        = cfg[3] ? :auto : cfg[4]        
	Emax        = cfg[5] ? nothing : cfg[6] 
	min_height  = cfg[7]
	min_prom    = cfg[8]
	smoothwin   = cfg[9]
end

# ╔═╡ 8d114fcf-17d6-43df-b2ed-7b2f01e918b2
begin

	records = DDR.load_eigenph_records(joinpath(data_root, mode))

	# ----- sort by Q -----
	sort!(records, by = r -> (ismissing(r.meta.Q) ? Inf : r.meta.Q))
	
	# -- preprocess; fix jumps in major resonances, and try to center
	#    the resonances so theyre all more or less in the same branch
	DDR.unwrap!(records; Emin=Emin === :auto ? nothing : Emin, period=pi)
	DDR.center_records_at!(records; Eref=0.1, period=pi)

	selector = selector_re

	# ----- resonance detection (use Emin_used/Emax_used) -----
	idxres = DDR.find_resonances_across_records(
		records;
	    selector=selector,
	    Emin=Emin,
	    Emax=Emax,
	    smoothwin=smoothwin,
	    min_height=min_height,
	    min_prominence=min_prom,
	)

	# ----- plot phases/derivatives with overlay -----
	plt = DDR.plot_eigenph(
		records;
	    selector=selector,
	    ytype=ytype,
	    smoothwin=smoothwin,
	    idxres=idxres,
		xlabel = "Energy (eV)",
		ylims= (ytype === :deriv || ytype === :abderiv) ? (0,200) : (:auto, :auto),

	)
	
	# ----- add the two vertical lines -----
	Emin === :auto   || vline!(plt, [Emin]; label=false, linestyle=:dash, alpha=0.6)
	Emax === nothing || vline!(plt, [Emax]; label=false, linestyle=:dash, alpha=0.6)
	(ytype === :absderiv || ytype === :deriv) && 
		hline!(plt, [min_height]; label=false, linestyle=:dash, alpha=0.5)
	
	plt
end

# ╔═╡ Cell order:
# ╠═9527d1fc-1f12-11f1-bbc3-ad70a776bfa7
# ╟─c47fdc16-11f3-4ac0-af04-f6468d48e989
# ╟─f6999a5f-9003-4270-ac1c-29a2350711ff
# ╟─94c35eb8-7e95-4f76-9829-94da62bc8878
# ╟─cc01839f-ed8a-46e4-b112-8d91404ea61a
# ╟─91911778-83ba-41c1-86a3-18ab1cb7c295
# ╟─a0854639-7af8-42e9-95f4-9561b7a6e420
# ╟─8d114fcf-17d6-43df-b2ed-7b2f01e918b2
