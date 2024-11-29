# the app.jl file is the main entry point to the application. It is a bridge between the UI and the data processing logic.
using GenieFramework, JLD2, DataFrames, CSV, StatsBase, PlotlyBase, DSP, ColorSchemes, Healpix, Random, Colors
@genietools

usvs_color_scheme = ColorSchemes.darkrainbow
max_pixels = 60

# const _model_filename = "20240912T190943289"
const _model_filename = "20241124T150609280"
const _available_eqs_jld2 = filter(x -> occursin(".jld2", x), readdir(joinpath("data", _model_filename), join=true))

# function load_eq_lat_long(available_eqs_jld2)
#     latv = Vector{Float32}()
#     lonv = Vector{Float32}()
#     depthv = Vector{Float32}()
#     for eqfile in available_eqs_jld2
#         eqfile = replace(eqfile, "$(_model_filename)" => "lat_lon_dep_data", ".jld2" => "_data.jld2")
#         eqfile = JLD2.jldopen(eqfile)
#         @show keys(eqfile)
#         lat = eqfile["latitude"]
#         lon = eqfile["longitude"]
#         depth = eqfile["depth"]
#         JLD2.close(eqfile)
#         push!(latv, lat)
#         push!(lonv, lon)
#         push!(depthv, depth)
#     end
#     return latv, lonv, depthv
# end

# const _eq_lat, _eq_long, _eq_depth = load_eq_lat_long(_available_eqs_jld2)

function find_closest_eq(lat, lon, eq_loc_data)
    distances = [sqrt((lat - eq_loc["latitude"])^2 + (lon - eq_loc["longitude"])^2) for eq_loc in eq_loc_data]
    return argmin(distances)
end

function jld2_to_eqname(jld2file)
    eqname = first((split(basename(jld2file), ".")))
    return eqname
end

envelope(x) = abs.(hilbert(x))

# const first_available_eq_name = available_eqs_name[1]
# default_eq = available_eqs_name[1]

# the measurements taken by each station are stored in a table
# df = DataFrame(Arrow.Table("data/data.lz4"))
# columns = names(df)
# the soil_data table contains the location of each station along with other information
# soil_data = CSV.read("data/soil_data.csv", DataFrame)
tgrid = collect(1:256)
stf = range(0.0, stop=1.0, length=256)

function read_all_JLD2_files(available_eqs_jld2)
    eqs = []
    for eqfile in available_eqs_jld2
        eqdat = JLD2.jldopen(eqfile)
        push!(eqs, eqdat)
    end
    return eqs
end

function read_all_loc_JLD2_files(available_eqs_jld2)
    eqs = []
    for eqfile in available_eqs_jld2
        eqfile = replace(eqfile, "$(_model_filename)" => "lat_lon_dep_data", ".jld2" => "_data.jld2")
        eqfile = JLD2.jldopen(eqfile)
        push!(eqs, eqfile)
    end
    return eqs
end

function get_min_KL_instance_id(jld2fileindex, eq_data)
    eqfile = eq_data[jld2fileindex]
    k = keys(eqfile)
    KLpixel = filter(x -> occursin("KL", x), collect(k))[1]
    id = argmin(eqfile[KLpixel][end, :])
    return id
end

function get_pixels(jld2fileindex, eq_data)
    eqfile = eq_data[jld2fileindex]
    k = keys(eqfile)
    pixels = filter(x -> tryparse(Int, x) !== nothing, collect(k))
    return pixels
end
function get_stf_bundle(jld2fileindex, eq_data, selected_pixels)
    eqfile = eq_data[jld2fileindex]
    stf = eqfile["$(selected_pixels)"]
    return stf
end

const _available_eqs = jld2_to_eqname.(_available_eqs_jld2)
const _eq_data = read_all_JLD2_files(_available_eqs_jld2)
const _eq_loc_data = read_all_loc_JLD2_files(_available_eqs_jld2)
@app begin
    @out available_eqs = _available_eqs
    @out eq_lat_selected = 0.0
    @out eq_long_selected = 0.0
    @in selected_eq = string(sample(_available_eqs))
    @out available_pixels = ["",]
    @onchange selected_eq begin
        jld_file_index = findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]
        available_pixels = get_pixels(jld_file_index, _eq_data)
        eq_lat_selected = _eq_loc_data[findfirst(x -> x == selected_eq, available_eqs)]["latitude"]
        eq_long_selected = _eq_loc_data[findfirst(x -> x == selected_eq, available_eqs)]["longitude"]
    end

@out eq_location_traces = [PlotlyBase.scattermapbox(    
                    lat=[eqloc["latitude"]],
                    lon=[eqloc["longitude"]],
                    mode="markers",
                    size=10,
                    marker=attr(size=10, color=eqloc["depth"], linecolor="black", colorscale="Blues", cmin=300, cmax=700),
                    name="",
                    hovertext=eqname) for (eqname, eqloc) in zip(_available_eqs, _eq_loc_data)]


    @out tab_ids = ["tab_stf", "tab_usvs"]
    @out tab_labels = ["Single Pixel", "All Pixels"]
    @in selected_tab = "tab_stf"


    @out tgrid = range(-200, stop=200, length=400)[101:356]
    @in selected_pixel = ""
    @onchange available_pixels begin
        selected_pixel = first(available_pixels)
        # notify(__model__, "Selected Pixel $selected_pixel")
    end
    # @out traces = [scatter(x=collect(1:10), y=randn(10)), scatter(x=collect(1:10), y=randn(10))]
    @out traces = [scatter()]
    @onchange selected_eq, selected_pixel begin
        jld_file_index = findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]
        stf_bundle = get_stf_bundle(jld_file_index, _eq_data, selected_pixel)
        stf = vec(envelope(stf_bundle["USVS"][:, get_min_KL_instance_id(jld_file_index, _eq_data)]))
        stf_std = vec(std(stf_bundle["USVS"], dims=2))
        stf_upper = stf .+ stf_std
        stf_lower = stf .- stf_std

        raw_env_mean = vec(stf_bundle["RAW"])
        raw_env_upper = vec(stf_bundle["RAW_upper_bound"])
        raw_env_lower = vec(stf_bundle["RAW_lower_bound"])

        # traces[!] = [] # With the [!] suffix we reassign the array without triggering a UI update
        traces = [
            scatter(
                x=tgrid,
                y=stf,
                line=attr(color="rgb(0,0,100)"),
                name="Source Time Function Using Variational SymAE",
                mode="lines",
                legendgroup="group1",
            ),
            scatter(
                x=vcat(tgrid, reverse(tgrid)), # x, then x reversed
                y=vcat(stf_upper, reverse(stf_lower)), # upper, then lower reversed
                fill="toself",
                fillcolor="rgba(0,0,100,0.5)",
                line=attr(color="rgba(0,0,0,0)"),
                hoverinfo="skip",
                legendgroup="group1",
                showlegend=false
            ),
            scatter(x=tgrid,
                y=raw_env_mean,
                line=attr(color="rgb(100,0,0)"),
                legendgroup="group2",
                name="Raw Displacement Seismogram Envelope",
                mode="lines"),
            scatter(
                x=vcat(tgrid, reverse(tgrid)), # x, then x reversed
                y=vcat(raw_env_upper, reverse(raw_env_lower)), # upper, then lower reversed
                fill="toself",
                fillcolor="rgba(100,0,0,0.2)",
                line=attr(color="rgba(0,0,0,0)"),
                legendgroup="group2",
                hoverinfo="skip",
                showlegend=false
            )
        ]
        # @push traces 
    end

    @out stf_layout = PlotlyBase.Layout(
        title=attr(
            x=0.5,                   # Horizontal position (0 is left, 0.5 is center, 1 is right)
            y=0.97,                   # Vertical position (0 is bottom, 1 is top)
            font=attr(size=22)     # Customize font size if needed
        ),
        template="plotly_white",
        height=700,
        # width=900,
        xaxis=attr(
            title="Time (s)",
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=22),
            nticks=10,
            gridwidth=1,
            range=(-50, 80),
        ),
        yaxis=attr(
            title="Normalized Amplitude",
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=22),
            nticks=10,
            gridwidth=1,
            range=(0, 8),
        ),
        legend=attr(
            orientation="h",         # Horizontal legend
            y=-0.2,                  # Position below x-axis (adjust this value as needed)
            font=attr(size=22),
            x=0.5,                   # Center the legend horizontally
            xanchor="center",        # Anchor the legend to the center of the x position
            yanchor="top",           # Align to the top so that it moves downward from x-axis
            traceorder="grouped",    # Group legend items to display on multiple rows
        )
    )


    @in selected_rec_pixel = ""
    @onchange available_pixels begin
        selected_rec_pixel = first(available_pixels)
    end
    @out eq_receiver_lat = []
    @out eq_receiver_lon = []
    @out eq_receiver_names = []
    @onchange selected_rec_pixel, selected_eq begin
        eq_loc_pixel = _eq_loc_data[findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]]["$(selected_rec_pixel)"]
        eq_receiver_lat = [eq_loc_pixel[k][1] for k in keys(eq_loc_pixel)]
        eq_receiver_lon = [eq_loc_pixel[k][2] for k in keys(eq_loc_pixel)]
        eq_receiver_names = [k for k in keys(eq_loc_pixel)]
    end

    @out traces_usvs = [scatter()]
    @onchange selected_eq, available_pixels, tgrid begin

        jld_file_index = findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]
        all_pixels = get_pixels(jld_file_index, _eq_data)
        traces_usvs[!] = [] # With the [!] suffix we reassign the array without triggering a UI update
        av = 0
        av_raw = 0
        scale = 2
        sorted_pixels = sort(tryparse.(Int, all_pixels))
        npixels = length(sorted_pixels)
        for (ipixel, pixel) in enumerate(sorted_pixels)
            angles = broadcast(pix2angRing(Resolution(4), pixel)) do x
                floor(Int, rad2deg(x))
            end
            jld_file_index = findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]
            stf_bundle = get_stf_bundle(jld_file_index, _eq_data, string(pixel))
            pixel_stf = vec(mean(envelope(stf_bundle["USVS"]), dims=2))
            pixel_raw = vec(stf_bundle["RAW"])

            av = av .+ pixel_stf
            av_raw = av_raw .+ pixel_raw

            current_yoffset = (scale * (ipixel - 1))
            pixel_stf = pixel_stf ./ maximum(pixel_stf) .* 2.5
            pixel_raw = pixel_raw ./ maximum(pixel_raw) .* 2.5

            # @show get(usvs_color_scheme, ipixel/npixels), ipixel/npixels
            push!(traces_usvs,
                scatter(
                    x=tgrid,
                    y=pixel_stf .+ current_yoffset,
                    # line = attr(color = "white"),#hex(get(usvs_color_scheme, ipixel/npixels))),
                    fillcolor=hex(get(usvs_color_scheme, ipixel / npixels)),
                    line=attr(color="black"),#hex(get(usvs_color_scheme, ipixel/npixels))),
                    fill="tonexty",
                    mode="lines",
                    name=string(pixel, ": ", angles),
                    xaxis="x2", yaxis="y",
                    legendgroup="group$ipixel",
                ))

            push!(traces_usvs,
                scatter(
                    x=tgrid,
                    y=pixel_raw .+ current_yoffset,
                    # line = attr(color = "white"),#hex(get(usvs_color_scheme, ipixel/npixels))),
                    fillcolor=hex(get(usvs_color_scheme, ipixel / npixels)),
                    line=attr(color="black"),#hex(get(usvs_color_scheme, ipixel/npixels))),
                    fill="tonexty",
                    mode="lines",
                    showlegend=false,
                    name=string(pixel, ": ", angles),
                    xaxis="x", yaxis="y",
                    legendgroup="group$ipixel",
                ))
        end
        last_yoffset = (scale * (32))
        av .= av ./ maximum(av) .* 2.5
        av_raw .= av_raw ./ maximum(av_raw) .* 2.5
        push!(traces_usvs,
            scatter(
                x=tgrid,
                y=av .+ last_yoffset,
                line=attr(color="black"),
                mode="lines",
                xaxis="x2", yaxis="y",
                name="mean",
                legendgroup="group0",
            ))
        push!(traces_usvs,
            scatter(
                x=tgrid,
                y=av_raw .+ last_yoffset,
                line=attr(color="black"),
                mode="lines",
                showlegend=false,
                xaxis="x", yaxis="y",
                name="mean",
                legendgroup="group0",
            ))
        @push traces_usvs # Update the traces vector when all traces are generated
    end


    @out usvs_layout = PlotlyBase.Layout(
        template="plotly_white",
        height=1200,
        title="Raw Envelope Stacking Vs. SymAE Source Time Functions",
        yaxis_anchor="x",
        yaxis=attr(showgrid=false, showticklabels=false),
        xaxis=attr(
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=22),
            range=(-50, 60),
            domain=[0, 0.48],
            showgrid=true,
            title="Time (s)"),
        xaxis2=attr(
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=22),
            range=(-50, 60),
            showgrid=true,
            domain=[0.52, 1],
            title="Time (s)"),
    )

    @out mapcolor = []
    @in data_click = Dict{String,Any}()  # data from map click event
    # when clicking on a plot, the data_click variable will be populated with the event info
    @onchange data_click begin
        if haskey(data_click, "points")
            lat_click, lon_click = data_click["points"][1]["lat"], data_click["points"][1]["lon"]
            closest_index = find_closest_eq(lat_click, lon_click, _eq_loc_data)
            selected_eq = available_eqs[closest_index]
        end

    end
end

# enable event detection in charts. Moreover, you'll need to add the class `sync_data` to the `plotly`
# element in the HTML code
@mounted watchplots()

# declare a route at / that'll render the HTML
@page("/", "app.jl.html")

