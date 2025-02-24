# the app.jl file is the main entry point to the application. It is a bridge between the UI and the data processing logic.
using GenieFramework, JLD2, DataFrames, CSV, StatsBase, PlotlyBase, DSP, ColorSchemes, Healpix, Random, Colors
using Dates, Distances
@genietools

# slab reading
const _slab_folders = filter(x -> occursin("Slab2_", x), readdir(joinpath("data", "Slab2"), join=true))

const _slab_dep_files = map(_slab_folders) do fld
    first(filter(x -> (occursin(".xyz", x) & occursin("dep", x)), readdir(fld, join=true)))
end

slab_dfs = map(_slab_dep_files) do file
    df = CSV.read(file, header=0, DataFrame)
end;

# source long is -180 to 180
function convert_longitude(lon)
    return lon > 180 ? lon - 360 : lon
end

# source long is -180 to 180
function find_nearest_slab(lat, long, slab_dfs)
    min_distance = Inf
    nearest_slab = nothing
    nearest_slab_index = nothing

    for (idf, df) in enumerate(slab_dfs)
        slab_lat = df[:, 2]
        slab_long = convert_longitude.(df[:, 1])

        # Calculate the Haversine distance between the given point and each point in the slab dataset
        distances = haversine.(Ref((lat, long)), zip(slab_lat, slab_long), 6372.8)
        # Calculate the mean distance
        distance = mean(distances)

        # Update the nearest slab and point if a closer one is found
        if distance < min_distance
            min_distance = distance
            nearest_slab = df
            nearest_slab_index = idf
        end
    end

    return nearest_slab, nearest_slab_index
end



# read background events
_background_events = CSV.read("data/deep_events_after2002_400km+.txt", DataFrame; delim='|')

function filter_background_events(df::DataFrame, given_date::Date, lat, lon)
    # Ensure the date column is parsed as DateTime type
    df.event_date = DateTime.(df[!, " Time "])

    df.latitude = df[!, " Latitude "]
    df.longitude = df[!, " Longitude "]
    # Define the date range
    start_date = given_date - Year(100)
    end_date = given_date + Year(100)

    # Define the latitude and longitude range
    d = 500

    # Filter the DataFrame using Haversine distance
    filtered_df = filter(row -> row.event_date >= start_date && row.event_date <= end_date &&
                                    haversine((lat, lon), (row.latitude, row.longitude), 6372.8) <= d, df)

    return filtered_df
end

const _usvs_color_scheme = ColorSchemes.darkrainbow
const _max_pixels = 60

# const _model_filename = "20240912T190943289"
# const _model_filename = "20241122T122806149"
# const _model_filename = "20241124T150609280"
# const _model_filename = "20250105T181423701"
# const _model_filename = "Lowest_mse_EachPixelOptim"
const _model_filename = "Mixed_Single_Multiple_Pixel_BestInterpreted_24Feb2025"



const _available_eqs_jld2 = filter(x -> occursin(".jld2", x), readdir(joinpath("data", _model_filename), join=true))

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
# const _tgrid = range(-200, stop=200, length=400)[101:356]
const _tgrid = range(-100, 100, length=200)

# stf = range(0.0, stop=1.0, length=256)

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

function get_stf_traces(selected_eq, selected_pixel, stf_enevelope_toggle)
    jld_file_index = findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]
    stf_bundle = _eq_data[jld_file_index]["$(string(selected_pixel))"]
    stf = dropdims(mean(envelope(stf_bundle["USVS"][:, :]), dims=2), dims=2)
    envelopes = envelope(stf_bundle["USVS"][:, :])
    mean_envelope = dropdims(mean(envelopes, dims=2), dims=2)

    # Calculate the correlation coefficient for each STF with the mean envelope
    correlations = [cor(envelopes[:, i], mean_envelope) for i in 1:size(envelopes, 2)]
    # Find the index of the STF with the highest correlation
    best_index = argmax(correlations)
    best_stf = stf_bundle["USVS"][:, best_index]
    stf = stf_enevelope_toggle ? mean_envelope : best_stf

    stf_std = dropdims(std(envelope(stf_bundle["USVS"][:, :]), dims=2), dims=(2))
    stf_upper = stf .+ stf_std
    stf_lower = stf .- stf_std

    raw_env_mean = (stf_bundle["RAW"])
    raw_env_upper = (stf_bundle["RAW_upper_bound"])
    raw_env_lower = (stf_bundle["RAW_lower_bound"])

    nw = length(stf) ÷ 30
    fs = step(_tgrid)
    spec = spectrogram(best_stf, nw, nw ÷ 2; fs=fs)
    spec_power = spec.power ./ maximum(spec.power)
    spec_tgrid = range(minimum(_tgrid), stop=maximum(_tgrid), length=size(spec.power, 2))
    # c = wavelet(Morlet(π), averagingType=NoAve(), β=3);
    # cwt_res, scales = cwt(stf, c)
    # frequencies = scal2frq(scales, c, inv(step(_tgrid)))

    # traces[!] = [] # With the [!] suffix we reassign the array without triggering a UI update
    traces = [
        scatter(
            x=_tgrid,
            y=stf,
            line=attr(color="rgb(0,0,128)"),
            name="Source Time Function Using Variational SymAE",
            mode="lines",
            legendgroup="group1",
            xaxis="x", yaxis="y",
        ),
        scatter(
            x=vcat(_tgrid, reverse(_tgrid)), # x, then x reversed
            y=vcat(stf_upper, reverse(stf_lower)), # upper, then lower reversed
            fill="toself",
            fillcolor="rgba(0,0,128,0.2)",
            line=attr(color="rgba(0,0,0,0)"),
            hoverinfo="skip",
            legendgroup="group1",
            showlegend=false,
            xaxis="x", yaxis="y",
        )]
    if (stf_enevelope_toggle)
        push!(traces,
            scatter(x=_tgrid,
                y=raw_env_mean,
                line=attr(color="rgb(139,0,0)"),
                legendgroup="group2",
                name="Raw Displacement Seismogram Envelope",
                xaxis="x", yaxis="y",
                mode="lines"))

        push!(traces,
            scatter(
                x=vcat(_tgrid, reverse(_tgrid)), # x, then x reversed
                y=vcat(raw_env_upper, reverse(raw_env_lower)), # upper, then lower reversed
                fill="toself",
                fillcolor="rgba(139,0,0,0.2)",
                line=attr(color="rgba(0,0,0,0)"),
                legendgroup="group2",
                hoverinfo="skip",
                xaxis="x", yaxis="y",
                showlegend=false
            ))
    end
    push!(traces,
        heatmap(x=spec_tgrid, y=spec.freq, z=spec_power, xaxis="x", yaxis="y2", colorscale="Earth"))
    return traces
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
        marker=attr(size=10, color=[eqloc["depth"]], autocolorscale=false, colorscale="Reds", cmin=300, cmax=600),
        name="",
        hovertext=eqname) for (eqname, eqloc) in zip(_available_eqs, _eq_loc_data)]


    @out tab_ids = ["tab_stf", "tab_usvs"]
    @out tab_ids2 = ["tab_receivers", "tab_sources"]
    @out tab_labels2 = ["Receiver Locations", "Regional Catalog Events And Slab 2.0"]
    @in selected_tab2 = "tab_receivers"
    @out tab_labels = ["Single Pixel", "All Pixels"]
    @in selected_tab = "tab_stf"



    @in slab_toggle = false
    @onchange selected_eq begin
        slab_toggle = false
    end
    # Create a 3D scatter plot
    @out background_eq_traces = [PlotlyBase.scatter3d()]
    @onchange selected_eq, slab_toggle begin
        given_date = Date(first(split(selected_eq, "_")), "yyyymmdd")
        given_lat = _eq_loc_data[findfirst(x -> x == selected_eq, available_eqs)]["latitude"]
        given_lon = _eq_loc_data[findfirst(x -> x == selected_eq, available_eqs)]["longitude"]
        given_depth = _eq_loc_data[findfirst(x -> x == selected_eq, available_eqs)]["depth"]


        filtered_event_df = filter_background_events(_background_events, given_date, given_lat, given_lon)
        event_latitudes = filtered_event_df[!, " Latitude "]
        event_longitudes = filtered_event_df[!, " Longitude "]
        event_depths = filtered_event_df[!, " Depth/km "]

        # println(given_lat, given_lon, given_depth)

        background_eq_traces = [
            PlotlyBase.scatter3d(
                x=event_longitudes,
                y=event_latitudes,
                z=event_depths,
                mode="markers",
                name="Regional Catalog Events",
                marker=attr(size=5, color="gray", showscale=false)
            ),
            PlotlyBase.scatter3d(
                x=[given_lon],
                y=[given_lat],
                z=[given_depth],
                mode="markers",
                name="Selected Event",
                marker=attr(size=5, color="black")
            )
        ]
        if (slab_toggle)
            slab_df, islab_df = find_nearest_slab(given_lat, given_lon, slab_dfs)
            slab_name = join(split(basename(_slab_folders[islab_df]), "_")[3:end]) * " slab"

            background_eq_traces = vcat(background_eq_traces, [
                PlotlyBase.scatter3d(
                    x=convert_longitude.(slab_df[:, 1]),
                    y=slab_df[:, 2],
                    z=-1.0 .* slab_df[:, 3],
                    mode="markers",
                    name=slab_name,
                    marker=attr(size=2, opacity=0.5, color=-1.0 .* slab_df[:, 3], colorscale="Reds", showscale=true, cmin=0, cmax=700)
                )
            ])
        end
    end

    @out background_eq_layout = PlotlyBase.Layout(
        scene=attr(
            xaxis=attr(title="Longitude"),
            yaxis=attr(title="Latitude"),
            zaxis=attr(title="Depth", autorange="reversed")
        ),
        legend=attr(
            orientation="h",  # Horizontal legend
            x=0.5,  # Center the legend horizontally
            y=-0.2,  # Position the legend below the plot
            xanchor="center",
            yanchor="top"
        ),
        coloraxis=attr(
            colorbar=attr(
                x=-1.1,  # Position the colorbar to the right of the plot
                y=0.5,  # Center the colorbar vertically
                yanchor="middle"
            )
        )
    )



    @in selected_pixel = ""
    @onchange available_pixels begin
        selected_pixel = first(available_pixels)
        # notify(__model__, "Selected Pixel $selected_pixel")
    end
    @in stf_envelope_toggle = false
    # @out traces = [scatter(x=collect(1:10), y=randn(10)), scatter(x=collect(1:10), y=randn(10))]
    @out stf_layout = PlotlyBase.Layout(
        title=attr(
            x=0.5,                   # Horizontal position (0 is left, 0.5 is center, 1 is right)
            y=0.97,                   # Vertical position (0 is bottom, 1 is top)
            font=attr(size=22)     # Customize font size if needed
        ),
        template="plotly_white",
        height=900,
        # width=900,
        ticklabelposition="inside top",
        xaxis=attr(
            title="Relative Time (s)",
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=15),
            nticks=10,
            showgrid=true,
            mirror=true,
            gridwidth=1,
            gridcolor="black",
            automargin=true,
            range=(-50, 80),
        ),
        yaxis=attr(
            title="Normalized Amplitude",
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=15),
            nticks=10,
            showgrid=true,
            mirror=true,
            gridwidth=1,
            gridcolor="black",
            # ticklabelposition="inside left",
            automargin=true,
            domain=[0, 0.75],
        ),
        yaxis2=attr(
            title="Frequency",
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=15),
            nticks=10,
            showgrid=true,
            mirror=true,
            gridwidth=1,
            gridcolor="black",
            # ticklabelposition="inside left",
            domain=[0.75, 1],
            automargin=true,
            # overlaying="y",
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
    @out traces = [scatter()]
    @onchange selected_eq, selected_pixel, stf_envelope_toggle begin
        traces = get_stf_traces(selected_eq, selected_pixel, stf_envelope_toggle)
    end

    @out usvs_layout = PlotlyBase.Layout(
        template="plotly_white",
        height=1200,
        yaxis_anchor="x",
        legend=attr(title="pixel (colatitude, longitude)"),
        yaxis=attr(showgrid=false, showticklabels=false),
        xaxis=attr(
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=15),
            range=(-50, 60),
            domain=[0, 0.48],
            showgrid=true,
            mirror=true,
            gridwidth=1, gridcolor="black",
            title="Relative Time (s)"),
        xaxis2=attr(
            titlefont=attr(size=22),
            font=attr(size=22),
            tickfont=attr(size=15),
            range=(-50, 60),
            showgrid=true,
            mirror=true,
            domain=[0.52, 1],
            gridwidth=1, gridcolor="black",
            title="Relative Time (s)"),
        dragmode="drawopenpath",
        newshape_line_color="black",
        newshape_line_width=1,
        newshape_line_style="dash",
        modebar_add=["drawline",
            "drawopenpath",
            "drawclosedpath",
            "drawcircle",
            "drawrect",
            "eraseshape"
        ],
        annotations=[
            attr(text="Raw Displacement Envelope Stacking",
                xref="paper", yref="paper",
                xanchor="center",
                font=attr(size=22),
                x=0.25, y=1.03, showarrow=false),
            attr(text="SymAE Source Time Functions",
                xref="paper", yref="paper",
                xanchor="center",
                font=attr(size=22),
                x=0.75, y=1.03, showarrow=false)
        ]
    )


    @in selected_rec_pixel = []
    @onchange available_pixels begin
        selected_rec_pixel = available_pixels
    end
    @out traces_receiver = [PlotlyBase.scattermapbox()]
    @onchange selected_rec_pixel, selected_eq begin
        jld_file_index = findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]
        all_pixels = get_pixels(jld_file_index, _eq_data)
        sorted_pixels = sort(tryparse.(Int, all_pixels))
        ipixels = map(x -> findfirst(y -> y == x, string.(sorted_pixels)), selected_rec_pixel)
        eq_receiver_lat = []
        eq_receiver_lon = []
        eq_receiver_names = []
        eq_receiver_colors = []
        for (ipixel, pixel) in zip(ipixels, selected_rec_pixel)
            if (pixel !== nothing)
                eq_loc_pixel = _eq_loc_data[findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]]["$(pixel)"]
                push!(eq_receiver_lat, [eq_loc_pixel[k][1] for k in keys(eq_loc_pixel)]...)
                push!(eq_receiver_lon, [eq_loc_pixel[k][2] for k in keys(eq_loc_pixel)]...)
                push!(eq_receiver_names, ["pixel $pixel: " * k for k in keys(eq_loc_pixel)]...)
                push!(eq_receiver_colors, ["#" * hex(get(_usvs_color_scheme, ipixel / length(sorted_pixels))) for i in 1:length(keys(eq_loc_pixel))]...)
            end
        end
        traces_receiver = [PlotlyBase.scattermapbox(
            lat=eq_receiver_lat,
            lon=eq_receiver_lon,
            mode="markers",
            size=10,
            marker=attr(size=6, autocolorscale=false, cauto=false, color=eq_receiver_colors),
            name="",
            hovertext=eq_receiver_names)]
    end

    @out traces_usvs = [scatter()]
    @onchange selected_eq begin

        jld_file_index = findall(x -> occursin(selected_eq, x), _available_eqs_jld2)[1]
        all_pixels = get_pixels(jld_file_index, _eq_data)
        av = 0
        av_raw = 0
        scale = 2
        sorted_pixels = sort(tryparse.(Int, all_pixels))


        traces_usvs_all = []

        EQ = _eq_data[jld_file_index]
        pixel_stds = [mean(std(envelope(EQ[string(p)]["USVS"]), dims=2)) for p in sorted_pixels]
        good_pixels = findall(x -> x < 0.25, pixel_stds)
        sorted_pixels = sorted_pixels[good_pixels]
        npixels = length(sorted_pixels)
        for (ipixel, pixel) in enumerate(sorted_pixels)
            angles = broadcast(pix2angRing(Resolution(4), pixel)) do x
                floor(Int, rad2deg(x))
            end
            stf_bundle = EQ["$(string(pixel))"]

            pixel_stf = dropdims(mean(envelope(stf_bundle["USVS"]), dims=2), dims=2)
            # pixel_stf = vec(mean(envelope(stf_bundle["USVS"]), dims=2))
            pixel_raw = stf_bundle["RAW"]

            av = av .+ pixel_stf
            av_raw = av_raw .+ pixel_raw

            current_yoffset = (scale * (ipixel - 1))
            pixel_stf = pixel_stf ./ maximum(pixel_stf) .* 2.5
            pixel_raw = pixel_raw ./ maximum(pixel_raw) .* 2.5
            # @show get(_usvs_color_scheme, ipixel/npixels), ipixel/npixels
            push!(traces_usvs_all,
                scatter(
                    x=_tgrid,
                    y=pixel_stf .+ current_yoffset,
                    # line = attr(color = "white"),#hex(get(_usvs_color_scheme, ipixel/npixels))),
                    fillcolor=hex(get(_usvs_color_scheme, ipixel / npixels)),
                    line=attr(color="black"),#hex(get(_usvs_color_scheme, ipixel/npixels))),
                    fill="tonexty",
                    mode="lines",
                    name=string(pixel, ": ", angles),
                    xaxis="x2", yaxis="y",
                    legendgroup="group$ipixel",
                ))

            push!(traces_usvs_all,
                scatter(
                    x=_tgrid,
                    y=pixel_raw .+ current_yoffset,
                    # line = attr(color = "white"),#hex(get(_usvs_color_scheme, ipixel/npixels))),
                    fillcolor=hex(get(_usvs_color_scheme, ipixel / npixels)),
                    line=attr(color="black"),#hex(get(_usvs_color_scheme, ipixel/npixels))),
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
        push!(traces_usvs_all,
            scatter(
                x=_tgrid,
                y=av .+ last_yoffset,
                line=attr(color="black"),
                mode="lines",
                xaxis="x2", yaxis="y",
                name="mean",
                legendgroup="group0",
            ))
        push!(traces_usvs_all,
            scatter(
                x=_tgrid,
                y=av_raw .+ last_yoffset,
                line=attr(color="black"),
                mode="lines",
                showlegend=false,
                xaxis="x", yaxis="y",
                name="mean",
                legendgroup="group0",
            ))
        traces_usvs = traces_usvs_all
        @push traces_usvs # Update the traces vector when all traces are generated
    end




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
@page("/", "app.jl.html", layout = Stipple.ReactiveTools.DEFAULT_LAYOUT(title="EQ SymAE"))

