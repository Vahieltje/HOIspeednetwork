using CairoMakie

######################################################################
#                              PLOTTING                              #
######################################################################
function make_abundance_graph(n_timeseries; is_modifier=false, plot_range::AbstractRange, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)
    Amount_plot_points = Int(1e5)
    if !is_modifier
        n_timeseries = n_timeseries[:, plot_range]
        n_timeseries = n_timeseries'
    else
        n_timeseries = n_timeseries[plot_range]
    end
    if !is_modifier
        num_its, N = size(n_timeseries)
        println("N = $N, num_its = $num_its")
    else
        num_its = length(n_timeseries)
    end

    if num_its >= Amount_plot_points
        plot_step = div(num_its, Amount_plot_points)
    else
        plot_step = 1
    end
    plot_timeseries = n_timeseries[1:plot_step:end, :]
    title = "species_abundances"
    ylabel = "n"
    if is_modifier
        title = "Modifier"
        ylabel = "m"

    end
    fig = Figure()
    ax = Axis(fig[1, 1]; plot_kwargs...)
    if !is_modifier
        for n in 1:size(plot_timeseries)[2]
            lines!(plot_timeseries[:, n], label="species $n")
        end
    else
        lines!(plot_timeseries[:])
    end
    if legend
        Legend(fig[1, 2], ax)
    end
    hidexdecorations!(ax, label=false)
    ax.xlabel = "t"
    ax.ylabel = ylabel
    Makie.ylims!(ax, -0.1, min(maximum(plot_timeseries) * 1.1, 10^10))
    if is_modifier
        Makie.ylims!(ax, minimum(plot_timeseries) - 1, min(maximum(plot_timeseries) * 1.1, 10^10))
    end
    # CairoMakie.plot(x, y; title = title, ylabel = ylabel, xlabel = "time") # label = legend
    return fig
end

function make_abundance_and_modifier_graph(n_timeseries::Matrix{Float64}, modifier_timeseries::Vector{Float64}; plot_range::AbstractRange, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)
    n_timeseries = n_timeseries[:, plot_range]
    n_timeseries = n_timeseries'

    modifier_timeseries = modifier_timeseries[plot_range]
    # modifier_timeseries = modifier_timeseries'
    num_its, N = size(n_timeseries)
    println("N = $N, num_its = $num_its")

    if num_its >= amount_plot_points
        plot_step = div(num_its, amount_plot_points)
    else
        plot_step = 1
    end
    plot_abundances_timeseries = n_timeseries[1:plot_step:end, :]
    plot_modifier_timeseries = modifier_timeseries[1:plot_step:end]

    fig = Figure()
    grid = fig[:, 1] = GridLayout()
    ax1 = Axis(grid[1, 1]; plot_kwargs...)
    ax2 = Axis(grid[2, 1])
    for n in 1:size(plot_abundances_timeseries)[2]
        lines!(ax1, plot_abundances_timeseries[:, n], label="species $n")
    end
    lines!(ax2, plot_modifier_timeseries)
    rowgap!(grid, 5)
    hidexdecorations!(ax1)
    if legend
        Legend(fig[:, 2], ax1)
    end
    ax1.xlabel = "t"
    ax1.ylabel = "n"
    ax2.xlabel = "t"
    ax2.ylabel = "m"
    ymax_abundances = min(maximum(plot_abundances_timeseries) * 1.1, 10^10)
    if ymax_abundances == NaN
        ymax_abundances = 10^10
    end
    #Makie.ylims!(ax1, -0.1, ymax_abundances)
    Makie.ylims!(ax2, minimum(plot_modifier_timeseries) - 1, max(maximum(plot_modifier_timeseries) * 1.1, 1))
    # CairoMakie.plot(x, y; title = title, ylabel = ylabel, xlabel = "time") # label = legend
    return fig
end

function make_abundance_and_modifier_graph(sd::Simulation_data; plot_range::AbstractRange=1:sd.last_it, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)
    make_abundance_and_modifier_graph(sd.n_timeseries, sd.m_timeseries;
        plot_range=plot_range, amount_plot_points=amount_plot_points, legend=legend, plot_kwargs...)
end

#TO DO: implement this for full_HOI_Simulation_data
function make_abundance_and_modifier_graph(sd::full_HOI_Simulation_data, m_line::Tuple{Int64,Int64}; plot_range::AbstractRange=1:sd.last_it, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)
    make_abundance_and_modifier_graph(sd.n_timeseries, sd.m_timeseries[m_line[1], m_line[2], :];
        plot_range=plot_range, amount_plot_points=amount_plot_points, legend=legend, plot_kwargs...)
end

function make_abundance_and_modifier_graph(sd::invasion_Simulation_data{one_HOI_Simulation_data}; plot_range::AbstractRange=1:sd.saved_steps, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)
    n_timeseries = get_n_timeseries(sd)
    m_timeseries = get_m_timeseries(sd)
    make_abundance_and_modifier_graph(n_timeseries, m_timeseries;
        plot_range=plot_range, amount_plot_points=amount_plot_points, legend=legend, plot_kwargs...)
end

function make_abundance_and_modifier_graph(sd::invasion_Simulation_data{full_HOI_Simulation_data}, m_line::Tuple{Int,Int}; plot_range::AbstractRange=1:sd.saved_steps, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)
    n_timeseries = get_n_timeseries(sd)
    m_timeseries = get_m_timeseries(sd)[m_line[1], m_line[2], :]
    make_abundance_and_modifier_graph(n_timeseries, m_timeseries;
        plot_range=plot_range, amount_plot_points=amount_plot_points, legend=legend, plot_kwargs...)
end

function make_abundance_and_nm_graph(n_timeseries::Matrix{Float64}, modifier_timeseries::Vector{Float64}, ijk_species::Tuple{Int64,Int64,Int64};
    plot_range::AbstractRange, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)

    nm_timeseries = n_timeseries[ijk_species[2], :] .* modifier_timeseries
    make_abundance_and_modifier_graph(n_timeseries, nm_timeseries;
        plot_range=plot_range, amount_plot_points=amount_plot_points, legend=legend, plot_kwargs...)
end

function make_abundance_and_nm_graph(sd::Simulation_data; plot_range::AbstractRange=1:sd.last_it, amount_plot_points::Int64=10^5, legend=false, plot_kwargs...)
    make_abundance_and_nm_graph(sd.n_timeseries, sd.m_timeseries, sd.ecosystem.ijk_species;
        plot_range=plot_range, amount_plot_points=amount_plot_points, legend=legend, plot_kwargs...)
end
