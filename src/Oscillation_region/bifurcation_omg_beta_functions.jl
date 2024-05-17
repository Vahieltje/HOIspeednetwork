using DrWatson
@quickactivate "Nataha community model"

include(srcdir("parameter_generation.jl"))
include(srcdir("simulation_functions.jl"))
include(srcdir("ugent_colors.jl"))

using CairoMakie
# using JLD2
using CSV
# using Tables
using DataFrames
using Distributed

# using Alert

function make_bifurcation_matrix(ω_range, β_range, ijk_species; αs, HOI_type, intransitive, minimum_abundance=10^(-7), stop_if_convergence=true)
  max_time_unmodified = 5 * 10^4
  dt = 3 / 10^3

  convergence_matrix = zeros(length(ω_range), length(β_range))
  coexistence_matrix = zeros(length(ω_range), length(β_range))
  N = 3

  Threads.@threads for i in 1:length(ω_range)
    #status = pmap(1:length(ω_range)) do i
    # for i in 1:length(ω_range)
    ω = 10^ω_range[i]

    # adapt the max calculation time to the speed of the system
    max_time = max(Int(round(max_time_unmodified / ω)), max_time_unmodified)
    # max_time = 10^7
    saved_steps = 10
    println("ω = $ω")
    pmap(1:length(β_range)) do j
      β = β_range[j]
      # print("β = $β, ")
      n, m, A, R, death_rates = generate_nonuniform_3_cycle(αs; intransitive=intransitive)
      n .= 1
      m[ijk_species[1], ijk_species[2]] = 1
      ecosys = HOI_type(N, n, m, A, β, ijk_species, R, death_rates, ω)
      sd = Simulation_data(ecosys; saved_timesteps=saved_steps, max_time, dt)
      simulate_society!(sd; minimum_abundance, stop_if_convergence)
      final_abundances = sd.n_timeseries[:, sd.last_it]
      coexistence = count(x -> (x > minimum_abundance), final_abundances)
      convergence_matrix[i, j] = sd.convergence
      coexistence_matrix[i, j] = coexistence
    end
  end
  return convergence_matrix, coexistence_matrix
end

"""
Zoek singuliere punten in de oscillatie/convergentie matrix
"""
function find_singular_nonconvergence(m, ω_range, β_range)
  xs = size(m, 1)
  ys = size(m, 2)
  for (i, ω) in enumerate(ω_range)
    for (j, β) in enumerate(β_range)
      if 1 < i < xs && 1 < j < ys
        surrounded = sum(m[i-1:i+1, j-1:j+1])
        if m[i, j] == false && surrounded == 8
          println("non-convergence for ω = $ω, β = $β")
        end
      end
    end
  end
  return nothing
end

"""
Generate the name to save the heatmap datapoints
"""
function filepath(; data_type="oscillation", αs, HOI_type, ijk_species, intransitive)
  if !(data_type in ["oscillation", "coexistence"])
    error("data_type should be either 'oscillation' or 'coexistence'")
  end
  base_folder = "b2basics/omg_β_bifurcation/"
  transitivity_folder = "trans_cycle/"
  if intransitive
    transitivity_folder = "intrans_cycle/"
  end
  data_type *= "_"
  base_path = datadir(base_folder * transitivity_folder * data_type)
  full_path = base_path * filename_last_part(αs, HOI_type, ijk_species)
  if !isdir(base_path)
    mkpath(base_path)
  end
  return full_path
end


function filename_last_part(αs, HOI_type, ijk_species)
  endstring = "_as=$(αs)"
  if HOI_type == single_HOI_ecosystem
    endstring *= "_single_"
  else
    endstring *= "_double_"
  end
  for i in ijk_species
    endstring *= "$(i)"
  end
  return endstring * ".csv"
end