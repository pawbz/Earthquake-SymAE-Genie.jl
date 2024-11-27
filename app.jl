# the app.jl file is the main entry point to the application. It is a bridge between the UI and the data processing logic.
using GenieFramework, JLD2, DataFrames, CSV, StatsBase, PlotlyBase, DSP, ColorSchemes, Healpix, Random
@genietools


usvs_color_scheme = ColorSchemes.inferno
max_pixels = 60

const model_filename = "20240912T190943289"
available_eqs_jld2 = filter(x -> occursin(".jld2", x), readdir(joinpath("..", "scripts", "SavedResults", model_filename), join=true))
# model_filename = "para-2024-09-12T19:09:43.289-batchsize=16_beta=1.0_network_type=conv_nhidden_dense=1024_nlayer_dense=3_nsteps=100_nt=256_ntau=19_p=30_q=30_spatial_transformer_gamma=1000.0_transformer=spatial.jld2"

function generate_random_lat_long(n)
    latitudes = rand(-90.0:0.1:90.0, n)
    longitudes = rand(-180.0:0.1:180.0, n)
    return [(lat, lon) for (lat, lon) in zip(latitudes, longitudes)]
end

eq_lat_long = generate_random_lat_long(length(available_eqs_jld2))

function find_closest_eq(lat, lon, eq_lat_long)
    distances = [sqrt((lat - eq[1])^2 + (lon - eq[2])^2) for eq in eq_lat_long]
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


function get_pixels(jld2file)
    eqfile = JLD2.jldopen(jld2file)
    k = keys(eqfile)
    pixels = filter(x -> tryparse(Int, x) !== nothing, collect(k))
    JLD2.close(eqfile)
    return pixels
end
function get_stf_bundle(jld2file, selected_pixels)
    eqfile = JLD2.jldopen(jld2file)
    stf = eqfile["$(selected_pixels)"]
    JLD2.close(eqfile)
    return stf
end
# pixels = get_pixels(available_eqs_jld2[2])
# println(pixels)
# get_stf(available_eqs_jld2[2], pixels[1])
# eqfile = JLD2.jldopen(joinpath("data/saved_nuisance_optimizers_1/P_new", "optimal_nuisance_results_bin_max_receivers-0.1-20130524_145631_okt-P-2.jld2"))
# get_pixels(filter(x -> occursin(last(jld2_to_eqname.(available_eqs_jld2)), x), available_eqs_jld2)[1])
# in the reactive code block, we'll implement the logic to handle user interaction


const _available_eqs = jld2_to_eqname.(available_eqs_jld2)
const _eq_lat = [eq[1] for eq in eq_lat_long]
const _eq_long = [eq[2] for eq in eq_lat_long]
@app begin
    @out available_eqs = _available_eqs
    @out eq_lat = _eq_lat
    @out eq_long = _eq_long
    @in selected_eq = string(sample(_available_eqs))
    @out available_pixels = ["",]
    @onchange selected_eq begin
        # notify(__model__, "Selected Earthquake $selected_eq")
        available_pixels = get_pixels(filter(x -> occursin(selected_eq, x), available_eqs_jld2)[1])
    end

    @out tab_ids = ["tab_stf", "tab_usvs"]
    @out tab_labels = ["STF", "USVS"]
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
        stf_bundle = get_stf_bundle(filter(x -> occursin(selected_eq, x), available_eqs_jld2)[1], selected_pixel)
        stf = vec(mean(envelope(stf_bundle["USVS"]), dims=2))
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
                name="Extracted P-Source Time Function Using Variational SymAE",
                mode="lines"
            ),
            scatter(
                x=vcat(tgrid, reverse(tgrid)), # x, then x reversed
                y=vcat(stf_upper, reverse(stf_lower)), # upper, then lower reversed
                fill="toself",
                fillcolor="rgba(0,0,100,0.5)",
                line=attr(color="rgba(0,0,0,0)"),
                hoverinfo="skip",
                showlegend=false
            ),
            scatter(x=tgrid,
                y=raw_env_mean,
                line=attr(color="rgb(100,0,0)"),
                name="Raw Displacement Seismogram P Envelope",
                mode="lines"),
            scatter(
                x=vcat(tgrid, reverse(tgrid)), # x, then x reversed
                y=vcat(raw_env_upper, reverse(raw_env_lower)), # upper, then lower reversed
                fill="toself",
                fillcolor="rgba(100,0,0,0.2)",
                line=attr(color="rgba(0,0,0,0)"),
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
        template="plotly_dark",
        height=700,
        # width=900,
        xaxis=attr(
            title="Time (s)",
            titlefont=attr(size=22),
            tickfont=attr(size=22),
            nticks=10,
            gridwidth=1,
            font_color="black",
            range=(-50, 80),
        ),
        yaxis=attr(
            title=attr(text="Normalized Amplitude",
                # x=0.1,
                # y=0.5,
                font=attr(size=22)),
            tickfont=attr(size=22),
            nticks=10,
            gridwidth=1,
            font_color="black",
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

    @out traces_usvs = [scatter(x=collect(1:10), y=randn(10)), scatter(x=collect(1:10), y=randn(10))]
    @onchange selected_eq, tgrid begin
        all_pixels = get_pixels(filter(x -> occursin(selected_eq, x), available_eqs_jld2)[1])
        traces_usvs[!] = [] # With the [!] suffix we reassign the array without triggering a UI update
        av = 0
        scale = 2
        sorted_pixels = sort(tryparse.(Int, all_pixels))
        npixels = length(sorted_pixels)
        # color_samples = [get(usvs_color_scheme, i/(npixels)) for i in 0:npixels-1]
        for (ipixel, pixel) in enumerate(sorted_pixels)
            angles = broadcast(pix2angRing(Resolution(4), pixel)) do x
                floor(Int, rad2deg(x))
            end
            stf_bundle = get_stf_bundle(filter(x -> occursin(selected_eq, x), available_eqs_jld2)[1], string(pixel))
            pixel_stf = vec(mean(envelope(stf_bundle["USVS"]), dims=2))
            av = av .+ pixel_stf
            current_yoffset = (scale * (ipixel - 1))
            pixel_stf = pixel_stf ./ maximum(pixel_stf) .* 2.5
            # @show get(usvs_color_scheme, ipixel/npixels), ipixel/npixels
            push!(traces_usvs,
            scatter(
                x=tgrid,
                y=pixel_stf .+ current_yoffset,
        	    line = attr(color = get(usvs_color_scheme, ipixel/npixels)),
                fill="tonexty",
                mode="lines",
                name=string(pixel, ": ", angles), 
            ))     
        end
        last_yoffset = (scale * (30))
        av .= av ./ maximum(av) .* 2.5
        push!(traces_usvs,
        scatter(
            x=tgrid,
            y=av .+ last_yoffset,
            line = attr(color = "black"),
            mode="lines",
            name="mean",
        ))
        @push traces_usvs   # Update the traces vector when all traces are generated
    end


    @out usvs_layout = PlotlyBase.Layout(
        showlegend=true,
        # template="plotly_dark",
        font_family="Computer Modern",
        font_color="black",
        font=attr(family="Computer Modern", size=22),
        width = 600,
        height = 1200,
        yaxis_showgrid=false,
        xaxis_range=(-50, 60),
        xaxis=attr(
            title="Time (s)",
            standoff=0,  # Distance from the axis
            tickfont=attr(size=22),
            nticks=10,
            gridwidth=1,
            gridcolor="Gray",
            font_family="Computer Modern",
            font_color="black",
            # xaxis_range = (-50, 100),
        ),
        title=attr(
            font=attr(size=30),
            #text = "$(eqname)",
            x=0.5,  # Center title horizontally
            xanchor="center",
            y=0.95,  # Position title closer to the top
            yanchor="top"
        ),
        # title=attr(text=eqname,font = attr(size = 30), x=-100,),
        yaxis=attr(showticklabels=false),
        legend=attr(font=attr(size=18)),
        margin=attr(b=150),
    )

    @in selected_rec_pixel = ""
    @onchange available_pixels begin
        selected_rec_pixel = first(available_pixels)
    end
    # @private available_pixels = []
    # @in selected_pixels = ""
    # @out stf = get_stf(first(available_eqs_jld2), first(get_pixels(first(available_eqs_jld2))))
    # @out tgrid = collect(1:256)


    # we first declare reactive variables to hold the state of the interactive UI components
    # for the select menu, we need a list of station codes, and a variable to hold the selected station code
    # @in selected_fips = 1001
    # same for the metrics
    # @out metrics = ["PS", "TS", "WS10M", "PRECTOT"]
    # @in selected_metrics = ["PS", "TS", "WS10M"]
    # we expose a dataframe containing the data for the line plot configured in Genie Builder's visual editor
    # @out fips_data = DataFrame("PS_date" => Date[], "PS" => Float64[], "TS_date" => Date[], "TS" => Float64[], "T2M_date" => Date[], "T2M" => Float64[], "WS10M_date" => Date[], "WS10M" => Float64[])
    # @in N = 400
    # map plot data, plot configured in GB
    @out mapcolor = []
    @in data_click = Dict{String,Any}()  # data from map click event
    # metric shown on map
    # @in map_metric = "TS"
    # date of data shown on map
    # @in date = "2000-01-21"
    # when selecting a new station, metric  or N update the plot data. The isready variable is a pre-defined boolean that is set 
    # to true when the page finishes loading
    # @onchange isready, selected_eq begin
    #     processing = true
    #     # fips_data[!] = DataFrame() # with the [!] suffix, the reactive variable changes in value but the update is not sent to the browser
    #     # for m in selected_metrics
    #     #     # the lttb function resamples the data in the array to a set of N points
    #     #     idx, resampled = lttb(df[df.fips .== selected_fips, m], N)
    #     #     fips_data[!, m*"_date"] = df[df.fips .== selected_fips, :date][idx]
    #     #     fips_data[!, m] = resampled
    #     # end
    #     # @push fips_data # push the updated data to the browser
    #     processing = false
    # end
    # update the map when picking a new date or metric.
    # @onchange isready, date, map_metric begin
    #     @show "updating map with data from $date"
    #     # mapcolor = df[df.date .== Date(date), map_metric]
    # end
    # @onchange selected_eq,isready begin
    #     available_pixels = get_pixels(filter(x->occursin(selected_eq, x), available_eqs_jld2)[1])
    #     notify(__model__, "Selected New Earthquake", :warning)
    # end

    # when clicking on a plot, the data_click variable will be populated with the event info
    @onchange data_click begin
        # each chart generates different event data. The map click event has a "points" field with the lat and lon
        if haskey(data_click, "points")
            lat_click, lon_click = data_click["points"][1]["lat"], data_click["points"][1]["lon"]
            #find closest fip
            # distances = sqrt.((soil_data.lat .- lat_click).^2 + (soil_data.lon .- lon_click).^2)
            closest_index = find_closest_eq(lat_click, lon_click, eq_lat_long)
            selected_eq = available_eqs[closest_index]
        # else
            # clicking on a line chart generates event data with x and y coordinates of the clicked point
            # date = string(Date(unix2datetime(data_click["cursor"]["x"] / 1000))) # the x point data is in Unix time
            # date_data = df[(df.date .== date) .& (df.fips .== selected_fips), metrics] |> Matrix |> vec
            # closest_index = argmin(abs.(date_data .- data_click["cursor"]["y"]))
            # map_metric = metrics[closest_index]
        end

    end
    # @onchange isready begin
    #     notify(__model__, "Run locally to work with the full dataset. See README on GitHub for instructions.", :warning)
    # end
end

# enable event detection in charts. Moreover, you'll need to add the class `sync_data` to the `plotly`
# element in the HTML code
@mounted watchplots()

# declare a route at / that'll render the HTML
@page("/", "app.jl.html")

