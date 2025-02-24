### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c12f89c7-a038-4b02-9616-81704b24aee6
using Distances

# ╔═╡ 6afc6078-702a-4cc1-8b12-4c1ab60902b8
using Statistics

# ╔═╡ 895fff1a-5c45-4e64-8e32-46fa57fad5f2
begin 
	
	using PlutoUI 
	
					using PlutoPlotly 
	using CSV
using DataFrames 
end

# ╔═╡ d6f7138c-2f6c-4023-9db8-f5a5313cc597
TableOfContents()

# ╔═╡ 9dda2f2e-959a-4b06-a403-914e38888bf7
md"# Read Directory"

# ╔═╡ 813a3390-5e31-490b-8689-6842495dfaf1
md" This is for calabria " 

# ╔═╡ 849a545b-d06b-4403-913b-e5ca587f8966
begin
	dir_path = "./data/Slab2"  # Directory
	folder = readdir(dir_path)
end

# ╔═╡ 52de43d1-e65c-47bb-a228-990157ca1d18
slab_folders = filter(x->occursin("Slab2_", x), readdir(joinpath("data", "Slab2"), join=true))

# ╔═╡ 15c3041a-d3f7-42ae-9f03-428184a189c4
@bind selected_folder confirm(MultiCheckBox(folder))

# ╔═╡ a01dbcc5-408e-4766-ba01-31fe87145b4c
slab_dep_files = map(slab_folders) do fld
	first(filter(x->(occursin(".xyz", x)&occursin("dep", x)), readdir(fld, join=true)))
	# joinpath("fld", )
end

# ╔═╡ 375a5244-9afd-48eb-80a3-7b489a4bb52f
begin
	dir_path_file="./data/Slab2/"*selected_folder[1]
     file=readdir(dir_path_file)
end

# ╔═╡ 1e4304a4-b11e-4d56-b06a-884f69ac5ea0
md" Here choose .xyz dep file. that is upper subducting plate, if you want thickness of slab then I have to use thk.xyz file also to add thickness at corresponding depth. "

# ╔═╡ c4b5c4af-7c86-4b8e-9f9a-cf11148b91c4
@bind selected_file confirm(MultiCheckBox(file))

# ╔═╡ a4882570-26b5-4bf3-9691-c9cc9c391591
slab_dfs = map(slab_dep_files) do file
	df=CSV.read(file,header=0, DataFrame)
end;

# ╔═╡ c7bfa6b9-b716-4256-aee0-136c35ae4bc9
slab_long = mapreduce(vcat, slab_dfs) do df
	df[:, 1]
end

# ╔═╡ 41e821b3-71d5-468b-b835-d18361803810
slab_lat = mapreduce(vcat, slab_dfs) do df
	df[:, 2]
end

# ╔═╡ a040324e-4b2a-4a8d-9c73-26af20c27a9b
# slab_depth = mapreduce(vcat, slab_dfs) do df
	# df[:, 3]
# end

# ╔═╡ 2f81a23e-405c-425d-8f44-b020ff7e386a
df=map(keys(selected_file)) do i 
	x=CSV.read("data/Slab2/"*selected_folder[i]*"/"*selected_file[i],header=0, DataFrame)
end

# ╔═╡ eaa202a6-3dce-4f67-8b3b-af886a9cf987
md" Read lat long and depth"

# ╔═╡ 62895cf7-fff2-49f4-a63a-5b37dd89b096
begin
	long=df[1][:,1]
	lat = df[1][:,2]
	depth=(df[1][:,3])
	
end

# ╔═╡ 65c0e73c-668e-4742-8004-7d54acf75bca
plot(depth)

# ╔═╡ b2e83063-5e24-4bb8-85e3-f6fba1f8f58a
findall(x->x>0, randn(10))

# ╔═╡ 158138b5-5c7b-4b29-a100-a1b6e296651c
md"# Plotting slab"

# ╔═╡ 6283c23a-ce92-11ef-16ee-23fb2172020b
begin
	scatter_plot = surface(
	    x = long,
	    y = lat,
		z=depth,
	    mode = "lines",
	    marker = attr(
	        # size = slip*0.1,  # Adjust the marker size
	        color = depth,  # Use time values to determine marker color
	        colorscale = "Hot",  # Choose a colormap
	        # colorbar = attr(title = "Time"),  # Add a colorbar with a title
			colorbar=true,
		),type=type="scatter3d"

	    # text = time,  # Display time values when hovering over points
	)
	
	# Create a layout with a map background
	layout = Layout(
	 scene=attr(xaxis_title="Lat",
        yaxis_title="Long",
        zaxis_title="Depth",
    
            xaxis=attr(autorange="reversed",
                nticks=3
               
            ),
            yaxis=attr(autorange="reversed",
                nticks=3
                
            ),
            zaxis=attr(
                nticks=3
                
            ),
	 ),
	    title = "Map with Time-Based Point Colors",
	)
	
	# Combine the plot and layout
	plot(scatter_plot, layout)
end

# ╔═╡ fdba8dd3-be75-4492-a1f6-f6e38b741170
begin
	
	a, b, d = 1.32, 1., 0.8
	c = a^2 - b^2
	dom = range(0, stop=2π, length=100)
	u = dom' .* ones(100)
	v = ones(100)' .* dom
	
	x = @. (d * (c - a * cos(u) * cos(v)) + b^2 * cos(u)) / (a - c * cos(u) * cos(v))
	y = @. b * sin(u) * (a - d*cos(v)) / (a - c * cos(u) * cos(v))
	z = @. b * sin(v) * (c*cos(u) - d) / (a - c * cos(u) * cos(v))
	
	p = make_subplots(
	    rows=1, cols=2,
	    specs=[Spec(kind="scene") Spec(kind="scene")],
	    subplot_titles=["Color corresponds to z"  "Color corresponds to distance to origin"]
	)
	add_trace!(p, surface(x=x, y=y, z=z, colorbar_x=-0.07), row=1, col=1)
	add_trace!(p, surface(x=long, y=lat, z=depth, surfacecolor=depth), row=1, col=2)
	relayout!(p, title_text="Ring cyclide")
	p
end

# ╔═╡ e9f111e9-978f-42ab-9a10-39e50bc01bc5
function nan_to_missing!(vec)
    return [isnan(v) ? missing : v for v in vec]
end


# ╔═╡ a408d2f4-eb16-4d1e-baf6-91f886bbc09a
plot(depth)

# ╔═╡ 14d94f6a-e0d9-4509-9d86-60d6e12dcab8
slab_dep_files[20]

# ╔═╡ dc6b1016-08ec-4dc0-9e4a-0bf850db7233
# source long is -180 to 180
function convert_longitude(lon)
    return lon > 180 ? lon - 360 : lon
end

# ╔═╡ ea4a0b07-a34d-41b7-939e-5ad994eee516

function find_nearest_slab(lat, long, slab_dfs)
    min_distance = Inf
    nearest_slab = nothing
    nearest_slab_index = nothing

    for (idf, df) in enumerate(slab_dfs)
        slab_lat = df[:, 2]
        slab_long = convert_longitude.(df[:, 1])
        
        # Calculate the Haversine distance between the given point and each point in the slab dataset
        distances = haversine.(Ref((lat, long)), zip(slab_lat, slab_long))
        # Calculate the mean distance
        distance = mean(distances)

        # Update the nearest slab and point if a closer one is found
        if distance < min_distance
			@show idf, mean(slab_lat), mean(slab_long)
            min_distance = distance
            nearest_slab = df
            nearest_slab_index = idf
        end
    end

    return nearest_slab, nearest_slab_index
end

# ╔═╡ cb511251-28ec-4b6f-bfe1-5959635e4363
find_nearest_slab(-13, -67, slab_dfs)

# ╔═╡ 62284599-cf73-4fa5-aa92-6d7f1d129782
let
	
	# Plotting slab
	heatmap_plot = densitymapbox(
	    lat = lat,
	    lon = long,
	    z = -1.0 * depth,
	    radius = 10,  # Adjust the radius for the heatmap
	    colorscale = "Reds",  # Choose a colormap
	    colorbar = attr(title = "Depth"),  # Add a colorbar with a title
	)
	
	# Create a layout with a map background
	layout = Layout(
	    mapbox = attr(
	        style = "open-street-map",
	        center = attr(lat = mean(lat), lon = mean(long)),
	        zoom = 5  # Adjust the zoom level
	    ),
	    title = "Heatmap of Slab Depths",
	)
	
	# Combine the plot and layout
	plot(heatmap_plot, layout)
end

# ╔═╡ cbe70971-41f9-4a90-a286-d61121a44a19
plot(randn(4))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.14"
DataFrames = "~1.6.1"
Distances = "~0.10.12"
PlutoPlotly = "~0.6.2"
PlutoUI = "~0.7.60"
Statistics = "~1.11.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "2b7257b558f1de2d8a9751f4689a193c647d8d54"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "545a177179195e442472a1c4dc86982aa7a1bef0"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.7"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "26ec26c98ae1453c692efded2b17e15125a5bea1"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.28.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

    [deps.Distances.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "7878ff7172a8e6beedd1dea14bd27c3c6340d361"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.22"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "9ebe25fc4703d4112cc418834d5e4c9a4b29087d"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.2"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "1147f140b4c8ddab224c94efa9569fc23d63ab44"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═d6f7138c-2f6c-4023-9db8-f5a5313cc597
# ╟─9dda2f2e-959a-4b06-a403-914e38888bf7
# ╟─813a3390-5e31-490b-8689-6842495dfaf1
# ╠═849a545b-d06b-4403-913b-e5ca587f8966
# ╠═52de43d1-e65c-47bb-a228-990157ca1d18
# ╠═15c3041a-d3f7-42ae-9f03-428184a189c4
# ╠═a01dbcc5-408e-4766-ba01-31fe87145b4c
# ╠═375a5244-9afd-48eb-80a3-7b489a4bb52f
# ╟─1e4304a4-b11e-4d56-b06a-884f69ac5ea0
# ╠═c4b5c4af-7c86-4b8e-9f9a-cf11148b91c4
# ╠═a4882570-26b5-4bf3-9691-c9cc9c391591
# ╠═c7bfa6b9-b716-4256-aee0-136c35ae4bc9
# ╠═41e821b3-71d5-468b-b835-d18361803810
# ╠═a040324e-4b2a-4a8d-9c73-26af20c27a9b
# ╠═2f81a23e-405c-425d-8f44-b020ff7e386a
# ╟─eaa202a6-3dce-4f67-8b3b-af886a9cf987
# ╠═62895cf7-fff2-49f4-a63a-5b37dd89b096
# ╠═65c0e73c-668e-4742-8004-7d54acf75bca
# ╠═b2e83063-5e24-4bb8-85e3-f6fba1f8f58a
# ╟─158138b5-5c7b-4b29-a100-a1b6e296651c
# ╠═6283c23a-ce92-11ef-16ee-23fb2172020b
# ╠═fdba8dd3-be75-4492-a1f6-f6e38b741170
# ╠═e9f111e9-978f-42ab-9a10-39e50bc01bc5
# ╠═a408d2f4-eb16-4d1e-baf6-91f886bbc09a
# ╠═cb511251-28ec-4b6f-bfe1-5959635e4363
# ╠═c12f89c7-a038-4b02-9616-81704b24aee6
# ╠═14d94f6a-e0d9-4509-9d86-60d6e12dcab8
# ╠═dc6b1016-08ec-4dc0-9e4a-0bf850db7233
# ╠═ea4a0b07-a34d-41b7-939e-5ad994eee516
# ╠═62284599-cf73-4fa5-aa92-6d7f1d129782
# ╠═cbe70971-41f9-4a90-a286-d61121a44a19
# ╠═6afc6078-702a-4cc1-8b12-4c1ab60902b8
# ╠═895fff1a-5c45-4e64-8e32-46fa57fad5f2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
