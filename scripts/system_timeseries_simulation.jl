using DrWatson
@quickactivate "HOIspeednetwork"

include(srcdir("simulation_functions.jl"))
include(srcdir("parameter_generation.jl"))
include(srcdir("plot_functions.jl"))

println("start")


##############
# parameters #
##############
#Ecological Parameters
N = 3                               # amount of species in the system
ω = 10^(0)                          # modifier speed
αs = [1,1,1]                        # pairwise interaction strengths: [α_AB, α_AC, α_BC]
β = -1.7                              # modification strength
HOI_species = (2, 1, 3)             # higher order interaction species: (i, j, k). i and j are the species involved in the pairwise interaction, k is the modifier
                                    # species 1 = A, 2= B, 3 = C
HOI_type = asymmetric_HOI_ecosystem     # higher order interaction type: single = asymmetric, double = symmetric
intransitive = true                 # is the system intransitive?


#Simulation settings
minimum_abundance = 1e-7            # minimum abundance to be considered present. If abundance goes lower, species is considered extinct
stop_if_convergence = false         # stop the simulation if the system converges to a stable state
max_time = 3*10^4                     # maximum number of timesteps of the simulation
dt = 3 / 10^3                       # timestep of the simulation       
saved_steps = 10^4                  # amount of timesteps to save in the timeseries. Lower = faster, higher = more accurate visualisation

# Plotting settings
logscale = false                  # plot the y-axis in logscale
savefig = true



##############
# Simulation #
##############
n, m, A, R, death_rates = generate_nonuniform_3_cycle(αs, intransitive=intransitive)

n.= [1, 1, 1]                            # initial abundances
m[HOI_species[1], HOI_species[2]] = 1    # initial modifier, only for the modified interaction
ecosys = HOI_type(N, n, m, A, β, HOI_species, R, death_rates, ω)    # create the ecosystem with the appropriate type of HOI

sd = Simulation_data(ecosys; saved_timesteps=saved_steps, max_time, dt) # create the simulation data object
sd.n_timeseries .= 0
println("start simulation")
@time simulate_society!(sd, 
    stop_if_convergence=stop_if_convergence,
    minimum_abundance=minimum_abundance)      # run the simulation
last_it, convergence = sd.last_it, sd.convergence   # get the number of iterations and the convergence status of the system at the end of the simulation


#################################
# Printing numbers for analysis #
#################################
final_abundances = sd.n_timeseries[:, last_it]
println("final abundances  = $final_abundances")
final_modifier = sd.m_timeseries[last_it]
println("final modifier = $final_modifier")
coexistence = count(x -> (x > minimum_abundance), final_abundances) / N
println("φ = $coexistence")

############
# Plotting #
############
presentation_theme = Theme(fontsize=30, linewidth=5)
set_theme!(presentation_theme)

plot_kwargs = ntuple2dict((title="αs = $αs, ω = $ω, β = $β, HOI_species = $HOI_species", xlabel="t", ylabel = "n"))

plot_kwargs = []

if logscale
  plot_kwargs = ntuple2dict((yscale=log10, ylabel="n", limits=(nothing, (1e-10,1.5))))
end

plot_range = 1:last_it
fig = make_abundance_and_modifier_graph(sd; plot_range=plot_range, amount_plot_points=saved_steps, legend=true, plot_kwargs...)
println("done")
ax = fig.content[1]
ylims!(ax, (1e-40, 1.7))


println("sd.convergence = $(sd.convergence)")


# Saves the plot to the plots folder
if savefig
  save_folder = "paper/"
  plot_name = "temp_timeseries.png"
  if !isdir(plotsdir(plot_name))
    mkpath(plotsdir(plot_name))
  end
  save(plotsdir(plot_name), fig)
end

fig
