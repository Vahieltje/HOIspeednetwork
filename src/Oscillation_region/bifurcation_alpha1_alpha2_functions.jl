using DrWatson
@quickactivate "Nataha community model"

include(srcdir("parameter_generation.jl"))
include(srcdir("simulation_functions.jl"))
include(srcdir("ugent_colors.jl"))

using CairoMakie
using JLD2
using CSV
# using Tables
using DataFrames

"""
Generates the matrix of oscillation and coexistence for the 3-cycle
"""
function make_bifurcation_matrix(ω_range, β_range, ijk_species; αs, HOI_type, intransitive)
  max_time = 10^7
  dt = 3 / 10^3
  saved_steps = 10^5
  minimum_abundance = 10^-7
  convergence_matrix = zeros(length(ω_range), length(β_range))
  coexistence_matrix = zeros(length(ω_range), length(β_range))
  N = 3

  for i in 1:length(ω_range)
    # for i in 1:length(ω_range)
    ω = 10^ω_range[i]
    for (j, β) in enumerate(β_range)
      n, m, A, R, death_rates = generate_nonuniform_3_cycle(αs; intransitive=intransitive)
      n .= 1
      m[ijk_species[1], ijk_species[2]] = 1
      ecosys = HOI_type(N, n, m, A, β, ijk_species, R, death_rates, ω)
      sd = Simulation_data(ecosys; saved_timesteps=saved_steps, max_time, dt)
      simulate_society!(sd; minimum_abundance)
      final_abundances = sd.n_timeseries[:, sd.last_it]
      coexistence = count(x -> (x > minimum_abundance), final_abundances)
      convergence_matrix[i, j] = sd.convergence
      coexistence_matrix[i, j] = coexistence
    end
  end
  return convergence_matrix, coexistence_matrix
end


"""
Generate the name to save the heatmap datapoints
"""
function filepath_αs(; data_type="oscillation", αs, HOI_type, ijk_species, intransitive)
  if !(data_type in ["oscillation", "coexistence"])
    error("data_type should be either 'oscillation' or 'coexistence'")
  end
  base_folder = "b2basics/"
  transitivity_folder = "trans_cycle/"
  if intransitive
    transitivity_folder = "intrans_cycle/"
  end
  data_type *= "_"
  base_path = datadir(base_folder * transitivity_folder * data_type)
  return base_path * filename_last_part(αs, HOI_type, ijk_species)
end


function filename_last_part(αs, HOI_type, ijk_species)
  endstring = "_as=$(αs)"
  if HOI_type == asymmetric_HOI_ecosystem
    endstring *= "_single_"
  else
    endstring *= "_double_"
  end
  for i in ijk_species
    endstring *= "$(i)"
  end
  return endstring * ".csv"
end

function filename_α(HOI_type, ijk_species, boosted_α="α12", intransitive=false)
  filename = "b2basics/"
  if intransitive
    filename *= "intrans_cycle/"
  else
    filename *= "trans_cycle/"
  end
  filename *= "ijk=$(ijk_species)_$(HOI_type)_$(boosted_α)_oscillations.csv"
  return datadir(filename)
end