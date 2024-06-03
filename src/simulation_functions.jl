using DrWatson
@quickactivate "Nataha community model"

using Pkg
using LoopVectorization
include(srcdir("simulation_datatypes.jl"))


"""
Contains the functions to run a simulation of a community, using only one HOI β for 3rd order interactions.
Uses the structs defined in simulation_datatypes.jl:
        HOI ecosystems to describe the community
        HOI simulation data to store the simulation data (contains the ecosystem)

"""
##########################################################################
#                   Ecological system simulation Main                    #
##########################################################################
#region 

"""
simulates a society of species. Does so in place.
Args:
    sd::Simulation_data: the simulation data containing the ecosystem to simulate
    stop_if_convergence::Bool: if true, the simulation will stop if the system converges to a stable state
    minimum_change_threshold::Float64: the minimum change in abundance to consider the system stable
    minimum_abundance::Float64: the minimum abundance to consider a species present

Returns:
    Nothing
"""
function simulate_society!(sd::Simulation_data; stop_if_convergence=true, minimum_change_threshold::Float64=10^-3 * sd.dt, minimum_abundance=10^-7)
    simulation_error_throwing(sd)

    max_time = sd.max_time
    dt = sd.dt
    n_timeseries = sd.n_timeseries
    m_timeseries = sd.m_timeseries
    n = sd.ecosystem.n
    m = sd.ecosystem.m
    N = length(n)

    maximum_change_threshold::Float64 = 10^10            # if the change in abundance is larger than this, the simulation will stop (exponential growth)

    # number to keep track to only save every few steps
    if max_time > sd.saved_timesteps
        save_interval = div(max_time, sd.saved_timesteps)
    else
        save_interval = 1
    end

    abs_abundance_change = ones(Float64, N)
    last_it = sd.saved_timesteps

    for timestep in 1:max_time
        update_state!(sd.ecosystem; dt=sd.dt, minimum_abundance=minimum_abundance)         # progress one timestep
        if timestep % save_interval == 0
            current_savepoint = div(timestep, save_interval)
            update_timeseries!(sd, current_savepoint)
            # Is a boolean to check if the simulation should be terminated.
            # Checks for convergence, and if the abundances are exploding.
            termination = termination_criteria!(abs_abundance_change, sd, current_savepoint; stop_if_convergence=stop_if_convergence, minimum_change_threshold=minimum_change_threshold, maximum_change_threshold=maximum_change_threshold)
            if termination == true
                last_it = current_savepoint
                break
            end
        end
    end
    finish_simulation!(sd, last_it)
end

#endregion
######################################################################
#                            Update state                            #
######################################################################
#region
function update_state!(h::HOI_ecosystem; dt::Real=1 / mean(h.R) / 1000, minimum_abundance=10^-4)
    n = h.n
    N = length(n)
    # Calculate dn/dt
    update_n!(h, dt)
    update_m!(h, dt)
    for i in 1:N
        if h.n[i] < minimum_abundance
            h.n[i] = 0.0
        end
    end
end

"""
updates n to the next timestep
n + n*(R - n + A.*m_2*n)*dt
"""
function update_n!(h::HOI_ecosystem, dt::Real=1 / mean(h.R) / 10000)
    N = length(h.n)
    n = h.n
    m = h.m
    A = h.A
    R = h.R
    d = h.d
    @avx for i = 1:N
        dni = R[i] + d[i] * n[i]
        for j = 1:N
            dni += A[i, j] * m[i, j] * n[j]
        end
        n[i] = n[i] * (1 + dni * dt)
    end
    return nothing
end

"""
updates m to the next timestep
mᵢⱼ + (βᵢⱼₖ*nₖ-mᵢⱼ-1)*dt
"""
function update_m!(ecosys::symmetric_HOI_ecosystem, dt::Float64)
    i = ecosys.HOI_species[1]
    j = ecosys.HOI_species[2]
    k = ecosys.HOI_species[3]
    dm = ecosys.β * ecosys.n[k] - ecosys.m[i, j] + 1
    ecosys.m[i, j] += dm * dt * ecosys.ω
    # makes sure m is symmetric
    ecosys.m[j, i] = ecosys.m[i, j]
    return nothing
end

function update_m!(h::asymmetric_HOI_ecosystem, dt::Float64)
    i = h.HOI_species[1]
    j = h.HOI_species[2]
    k = h.HOI_species[3]
    mᵢⱼ = h.m[i, j]
    dm::Float64 = h.n[k] * h.β - mᵢⱼ + 1
    Δm::Float64 = dm * dt * h.ω
    h.m[i, j] += Δm
    return nothing
end

function update_m!(h::full_HOI_ecosystem, dt::Float64)
    N = length(h.n)
    @avx for i = 1:N
        for j = 1:N
            dmij = 1 - h.m[i, j]
            for k = 1:N
                dmij += h.B[i, j, k] * h.n[k]
            end
            h.m[i, j] += dmij * dt * h.Ω[i, j]
        end
    end
    return nothing
end

#endregion
######################################################################
#                  Update timeseries & termination                   #
######################################################################
#region

function update_timeseries!(sd::one_HOI_Simulation_data, current_savepoint::Int)
    n_timeseries = sd.n_timeseries
    m_timeseries = sd.m_timeseries
    n = sd.ecosystem.n
    m = sd.ecosystem.m
    HOI_species = sd.ecosystem.HOI_species
    n_timeseries[:, current_savepoint] = n
    m_timeseries[current_savepoint] = m[HOI_species[1], HOI_species[2]]
    return nothing
end

function update_timeseries!(sd::full_HOI_Simulation_data, current_savepoint::Int)
    n_timeseries = sd.n_timeseries
    m_timeseries = sd.m_timeseries
    n = sd.ecosystem.n
    m = sd.ecosystem.m
    n_timeseries[:, current_savepoint] = n
    m_timeseries[:, :, current_savepoint] .= m
    return nothing
end

function termination_criteria!(abs_abundance_change::Vector{Float64}, sd::Simulation_data, current_savepoint::Int, ; minimum_change_threshold::Float64=10^-4, maximum_change_threshold::Float64=10^10, stop_if_convergence::Bool=true)
    if current_savepoint <= 1
        return false
    end
    n_timeseries = sd.n_timeseries
    n = sd.ecosystem.n
    @views abs_abundance_change .= abs.(n_timeseries[:, current_savepoint-1] .- n)
    @views max_abs_abundance_change = maximum(abs_abundance_change)
    # terminate if convergence
    if max_abs_abundance_change <= minimum_change_threshold
        sd.convergence = true
        if stop_if_convergence == true
            return true
        end
    end
    # terminate if abundance explodes
    if max_abs_abundance_change >= maximum_change_threshold
        return true
    end
    return false
end

function finish_simulation!(sd::one_HOI_Simulation_data, last_it)
    sd.n_timeseries[:, last_it] .= sd.ecosystem.n
    sd.n_timeseries[:, last_it+1:end] .= 0
    sd.m_timeseries[last_it] = sd.ecosystem.m[sd.ecosystem.HOI_species[1], sd.ecosystem.HOI_species[2]]
    sd.m_timeseries[last_it+1:end] .= 0
    sd.last_it = last_it
end


function finish_simulation!(sd::full_HOI_Simulation_data, last_it)
    """
    Sets the last timestep of the timeseries to the last timestep of the simulation, and sets everything thereafter to 0
    Necessary for repeated experiments, as the timeseries are reused.
    """
    n_timeseries = sd.n_timeseries
    m_timeseries = sd.m_timeseries
    n_timeseries[:, last_it] .= n
    n_timeseries[:, last_it+1:end] .= 0
    m_timeseries[:, :, last_it] .= m
    m_timeseries[:, :, last_it+1:end] .= 0
    sd.last_it = last_it
end

#endregion
######################################################################
#                        Invasion simulation                         #
######################################################################
#region
"""
Simulates an ecosystem with invasion
Args:
    sd::Simulation_data: the simulation data using an ecosystem where both the original species and invaders are present
        the invasive species is the last one in the ecoystem.n vector.
    invasion_point::Int: the time when invasion happens (after how many timesteps)
    invasion_amount::Float64: the amount of invaders added to the ecosystem
"""
function simulate_society!(invasion_sd::invasion_Simulation_data; stop_if_convergence=true)
    sd_pre_invasion = invasion_sd.sd_pre_invasion
    sd_post_invasion = invasion_sd.sd_post_invasion
    simulate_society!(sd_pre_invasion, stop_if_convergence=true)
    invasion_sd.sd_post_invasion.ecosystem.n[end] = invasion_sd.invasion_amount
    simulate_society!(sd_post_invasion; stop_if_convergence)
    invasion_sd.last_it = sd_pre_invasion.last_it + sd_post_invasion.last_it
    return nothing
end


function simulation_error_throwing(sd::one_HOI_Simulation_data)
    n_timeseries = sd.n_timeseries
    m_timeseries = sd.m_timeseries
    N = length(sd.ecosystem.n)
    if length(m_timeseries) != size(n_timeseries)[2]
        error("m2_timeseries should have same amount of timesteps as n_timeseries")
    end
    if size(n_timeseries)[1] != N
        error("n_timeseries should have dimension $N x timestep, but has dimensions $(size(n_timeseries))")
    end
    return nothing
end

function simulation_error_throwing(sd::full_HOI_Simulation_data)
    n_timeseries = sd.n_timeseries
    m_timeseries = sd.m_timeseries
    N = length(sd.ecosystem.n)
    if length(m_timeseries[1, 1, :]) != size(n_timeseries)[2]
        error("m2_timeseries should have same amount of timesteps as n_timeseries")
    end
    if size(n_timeseries)[1] != N
        error("n_timeseries should have dimension $N x timestep, but has dimensions $(size(n_timeseries))")
    end
    return nothing
end

#endregion
