"""
This file contains the structs to store the data of the ecosystems and the simulations.
"""

################################################################################
#                              Ecosystem structs                               #
################################################################################
abstract type HOI_ecosystem end
# a type for HOI_ecosystem that only allows one HOI interaction, either working on one or two ways of the interaction
abstract type one_HOI_ecosystem <: HOI_ecosystem end
"""
An ecosystem where all HOI interactions can be present.
"""
struct full_HOI_ecosystem <: HOI_ecosystem
    N::Int64
    n::Vector{Float64}
    m::Matrix{Float64}
    A::Matrix{Float64}
    B::Array{Float64,3}
    R::Vector{Float64}
    d::Vector{Float64}
    Ω::Matrix{Float64}
end

function full_HOI_ecosystem(; N::Int=20, µₐ::Real=-0.0, σₐ::Real=0.0, Γ::Real=0.0,
    µᵦ::Real=0.0, σᵦ::Real=0.0, µᵣ::Real=10^5, σᵣ::Real=0.0, µ_d::Real=-1, σ_d::Real=0.0, Ω::Real=10)
    return full_HOI_ecosystem(N, generate_parameters(N=N, µₐ=µₐ, σₐ=σₐ, Γ=Γ, µᵦ=µᵦ, σᵦ=σᵦ,
            µᵣ=µᵣ, σᵣ=σᵣ, µ_d=µ_d, σ_d=σ_d)..., Ω)
end

function full_HOI_ecosystem(h::full_HOI_ecosystem; Ω)
    return full_HOI_ecosystem(h.N, h.n, h.m, h.A, h.B, h.R, h.d, Ω)
end

function Base.copy(h::full_HOI_ecosystem)
    return full_HOI_ecosystem(h.N, copy(h.n), copy(h.m), copy(h.A), copy(h.B), copy(h.R), copy(h.d), h.Ω)
end


"""
An ecosystem where only one HOI interaction can be present. The interaction is specified by the tuple HOI_species.
"""
struct single_HOI_ecosystem <: one_HOI_ecosystem
    N::Int64
    n::Vector{Float64}
    m::Matrix{Float64}
    A::Matrix{Float64}
    β::Float64
    HOI_species::Tuple{Int64,Int64,Int64}
    R::Vector{Float64}
    d::Vector{Float64}
    ω::Float64
end

function Base.copy(h::single_HOI_ecosystem)
    return single_HOI_ecosystem(h.N, copy(h.n), copy(h.m), copy(h.A), copy(h.β), h.HOI_species, copy(h.R), copy(h.d), h.ω)
end

function single_HOI_ecosystem(h::single_HOI_ecosystem; ω)
    return single_HOI_ecosystem(h.N, h.n, h.m, h.A, h.β, h.HOI_species, h.R, h.d, ω)
end

"""
An ecosystem where only one HOI interaction can be present, but both directions of the interaction are influenced by the modifier. The interaction is specified by the tuple HOI_species.
"""
struct symmetric_HOI_ecosystem <: one_HOI_ecosystem
    N::Int64
    n::Vector{Float64}
    m::Matrix{Float64}
    A::Matrix{Float64}
    β::Float64
    HOI_species::Tuple{Int64,Int64,Int64}
    R::Vector{Float64}
    d::Vector{Float64}
    ω::Float64
end

function Base.copy(h::symmetric_HOI_ecosystem)
    return symmetric_HOI_ecosystem(h.N, copy(h.n), copy(h.m), copy(h.A), h.β, copy(h.HOI_species), copy(h.R), copy(h.d), h.ω)
end


################################################################################
#                           Simulation Data structs                            #
################################################################################

abstract type Simulation_data end


"""
A struct containing all the data pertaining to a simulation
ecosystem, abundances_timeseries, m_timeseries, last_it, convergence
"""
mutable struct one_HOI_Simulation_data{T<:one_HOI_ecosystem} <: Simulation_data
    const ecosystem::T
    const n_timeseries::Matrix{Float64}
    const m_timeseries::Vector{Float64}
    # Iteration where the simulation stopped because all change stopped, or divergence happened
    last_it::Int64
    # Whether or not the system converges (either fixed value or oscillations if true)
    convergence::Bool
    const saved_timesteps::Int64
    const max_time::Int64
    const dt::Float64
    function one_HOI_Simulation_data{T}(ecosystem; saved_timesteps=10^5, max_time=10^6, dt=1e-3) where {T<:one_HOI_ecosystem}
        n_timeseries = Array{Float64}(undef, ecosystem.N, saved_timesteps)
        m_timeseries = Vector{Float64}(undef, saved_timesteps)
        return new{T}(ecosystem, n_timeseries, m_timeseries, 0, false, saved_timesteps, max_time, dt)
    end
end


"""
A struct containing all the data pertaining to a simulation
ecosystem, abundances_timeseries, m_timeseries, last_it, convergence
"""
mutable struct full_HOI_Simulation_data <: Simulation_data
    const ecosystem::full_HOI_ecosystem
    n_timeseries::Matrix{Float64}
    m_timeseries::Array{Float64,3}
    # Iteration where the simulation stopped because all change stopped, or divergence happened
    last_it::Int
    # Whether or not the system converges (either fixed value or oscillations if true)
    convergence::Bool
    const saved_timesteps::Int64
    const max_time::Int64
    const dt::Float64
    function full_HOI_Simulation_data(ecosystem::full_HOI_ecosystem; saved_timesteps=10^5, max_time=10^6, dt=1e-3)
        abundances_timeseries = Array{Float64}(undef, N, saved_timesteps)
        m_timeseries = Array{Float64}(undef, N, N, saved_timesteps)
        return new(ecosystem, abundances_timeseries, m_timeseries, 0, false, saved_timesteps, max_time, dt)
    end
end

function Simulation_data(ecosystem::T; saved_timesteps=10^5, max_time=10^6, dt=1e-3) where {T<:one_HOI_ecosystem}
    return one_HOI_Simulation_data{T}(ecosystem; saved_timesteps=saved_timesteps, max_time=max_time, dt=dt)
end

function Simulation_data(ecosystem::T; saved_timesteps=10^5, max_time=10^6, dt=1e-3) where {T<:full_HOI_ecosystem}
    return full_HOI_Simulation_data(ecosystem; saved_timesteps=saved_timesteps, max_time=max_time, dt=dt)
end

#########################################################################
#                               Invasion                                #
#########################################################################

"""
A struct containing all the data pertaining to a simulation with invasion. Consists of 2 Simulation_data structs, one for the society before the invasion, and one for the society after the invasion.
ecosystem, abundances_timeseries, m_timeseries, last_it, convergence
"""
mutable struct invasion_Simulation_data{T<:Simulation_data}
    const sd_pre_invasion::T
    const sd_post_invasion::T
    const invasion_amount::Float64
    const saved_timesteps::Int64
    last_it::Int64
    function invasion_Simulation_data{T}(ecosystem, invasion_time_fraction, invasion_amount; saved_timesteps=10^5, max_time=10^6, dt=1e-3) where {T<:Simulation_data}
        saved_timesteps_pre = Int(round(saved_timesteps * invasion_time_fraction))
        saved_timesteps_post = saved_timesteps - saved_timesteps_pre
        max_time_pre = Int(round(max_time * invasion_time_fraction))
        max_time_post = max_time - max_time_pre
        sd_pre_invasion = Simulation_data(ecosystem; saved_timesteps=saved_timesteps_pre, max_time=max_time_pre, dt=dt)
        sd_post_invasion = Simulation_data(ecosystem; saved_timesteps=saved_timesteps_post, max_time=max_time_post, dt=dt)
        return new{T}(sd_pre_invasion, sd_post_invasion, invasion_amount, saved_timesteps, 0)
    end
end

function get_n_timeseries(sd::invasion_Simulation_data)
    @views n_timeseries = [sd.sd_pre_invasion.n_timeseries[:, 1:sd.sd_pre_invasion.last_it] sd.sd_post_invasion.n_timeseries[:, 1:sd.sd_post_invasion.last_it]]
    return n_timeseries
end

function get_m_timeseries(sd::invasion_Simulation_data{<:one_HOI_Simulation_data})
    @views m_timeseries = [sd.sd_pre_invasion.m_timeseries[1:sd.sd_pre_invasion.last_it]; sd.sd_post_invasion.m_timeseries[1:sd.sd_post_invasion.last_it]]
    return m_timeseries
end

function get_m_timeseries(sd::invasion_Simulation_data{<:full_HOI_Simulation_data})
    @views m_timeseries = cat(sd.sd_pre_invasion.m_timeseries[:, :, 1:sd.sd_pre_invasion.last_it], sd.sd_post_invasion.m_timeseries[:, :, 1:sd.sd_post_invasion.last_it]; dims=3)
    return m_timeseries
end