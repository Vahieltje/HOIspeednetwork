using DrWatson
@quickactivate "Nataha community model"

include(srcdir("Oscillation_region/bifurcation_omg_beta_functions.jl"))


"""
Generates the matrix and figure of oscillation and coexistence regions for the 3-species HOI systems 
"""


##############
# parameters #
##############

ω_range = -3:0.4:2
β_range = -80:4:0
# β_range = 0:1:40

α12 = 1
α23 = 1

αs = [1.0, 1.0, 1.0]
#αs = [1,0,1]

HOI_type = symmetric_HOI_ecosystem         # can be asymmetric_HOI_ecosystem or symmetric_HOI_ecosystem
ijk_species = (1, 2, 3)
intransitive = true
minimum_abundance = 1e-7
stop_if_convergence = true          # stops the simulation of one system if it converges to a stable state


#########################
# New data generation #
#########################

get_new_data = true             # if true, new data will be generated, if false, data will be loaded from file
save_new_data = false & get_new_data
savefig = true

########################
# getting heatmap data #
########################

println("αs = $αs, HOI_type = $HOI_type, ijk_species = $ijk_species, minimum_abundance = $minimum_abundance, intransitive = $intransitive")
oscillation_file = filepath(; data_type="oscillation", αs, HOI_type, ijk_species, intransitive)
coexistence_file = filepath(; data_type="coexistence", αs, HOI_type, ijk_species, intransitive)

if get_new_data
    convergence_matrix, coexistence_matrix = make_bifurcation_matrix(ω_range, β_range, ijk_species; αs, HOI_type, intransitive, minimum_abundance, stop_if_convergence)
else
    convergence_matrix = Matrix(CSV.read(oscillation_file, DataFrame))
    coexistence_matrix = Matrix(CSV.read(coexistence_file, DataFrame))
end


# ###########
# CSV save #
# ###########
if save_new_data
    CSV.write(oscillation_file, Tables.table(convergence_matrix), writeheader=false)
    CSV.write(coexistence_file, Tables.table(coexistence_matrix), writeheader=false)
end

############
# Plotting #
############

colormap_oscillation = [ugent_yellow, ugent_blue]
colormap_oscillation_inverted = [colormap_oscillation[2], colormap_oscillation[1]]
# Oscillation heatmap
println("plotting heatmap")
fig = Figure()
ax = Axis(fig[1, 1], xlabel="log ω", ylabel="β", aspect=1, title="α = $αs")
plot_range = 101:501
hm = Makie.heatmap!(ax, ω_range, β_range, convergence_matrix, colormap=colormap_oscillation, colorrange=(0, 1))

resize_to_layout!(fig)

# coexistence heatmap
coexistence_colormap = ["#000000", "#002200", "#244622", "#4f8925"]
fig2 = Figure()
ax2 = Axis(fig2[1, 1], xlabel="log ω", ylabel="β", aspect=1, title="αs = $αs")
hm2 = Makie.heatmap!(ax2, ω_range, β_range, coexistence_matrix, colorrange=(0, 3),
    colormap=cgrad(coexistence_colormap, 4, categorical=true))


# cb2 = Colorbar(fig2[1, 2], hm2, label="Surviving species")
# colsize!(fig2.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig2)


#######################################
# Overlay plot oscillation/coexitence #
#######################################
convergence_oscillation_matrix = coexistence_matrix .* convergence_matrix

mixed_colormap = [ugent_yellow, "#002200", "#244622", "#4f8925"]

fig3 = Figure()
ax3 = Axis(fig3[1, 1], xlabel="log ω", ylabel="β", aspect=1, title="αs = $αs")
hm2 = Makie.heatmap!(ax3, ω_range, β_range, convergence_oscillation_matrix, colorrange=(0, 3),
    colormap=cgrad(mixed_colormap, 4, categorical=true))



if savefig
    save_folder = "paper/"
    plot_name_1 = "osc_α=$(αs[1]).png"
    save(plotsdir(plot_name_1), fig)
    plot_name_2 = "coex_α=$(αs[1]).png"
    save(plotsdir(plot_name_2), fig2)
    plot_name_3 = "coex_osc_α=$(αs[1]).png"
    save(plotsdir(plot_name_3), fig3)
    println("saved $plot_name_3")
end

fig3
