using DrWatson
@quickactivate "Nataha community model"

include(srcdir("b2basics/bifurcation_omg_beta_functions.jl"))


"""
Generates the matrix and figure of oscillation and coexistence regions for the 3-species HOI systems 
"""


##############
# parameters #
##############

ω_range = -3:0.1:2
β_range = -80:1:0
# β_range = 0:1:40

α12 = 2
α23 = 2

αs = [α12, α23, α23]
#αs = [1,0,1]

HOI_type = single_HOI_ecosystem
ijk_species = (1, 2, 3)
intransitive = true
superiority = (1, 2, 3)
minimum_abundance = 1e-70
stop_if_convergence = false

get_new_data = true
save_new_data = false & get_new_data
savefig = true


########################
# getting heatmap data #
########################

println("αs = $αs, HOI_type = $HOI_type, ijk_species = $ijk_species, minimum_abundance = $minimum_abundance, intransitive = $intransitive")
oscillation_file = filepath(; data_type="oscillation", αs, HOI_type, ijk_species, intransitive, superiority)
coexistence_file = filepath(; data_type="coexistence", αs, HOI_type, ijk_species, intransitive, superiority)

if get_new_data
    convergence_matrix, coexistence_matrix = make_bifurcation_matrix(ω_range, β_range, ijk_species; αs, HOI_type, intransitive, superiority, minimum_abundance, stop_if_convergence)
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
ax = Axis(fig[1, 1], xlabel="modifier speed (log ω)", ylabel="Modification strength (β)", aspect=1, title="α = $αs")
plot_range = 101:501
hm = Makie.heatmap!(ax, ω_range, β_range, convergence_matrix, colormap=colormap_oscillation, colorrange=(0, 1))

# subgrid = GridLayout(fig[1, 2], tellheight=false)
# Label(subgrid[1, 1], "Oscillating")
# Label(subgrid[3, 1], "Not\nOscillating")
# cb = Colorbar(subgrid[2, 1], limits=(0, 1),
#     colormap=cgrad(colormap_oscillation_inverted, 2, categorical=true),
#     ticksvisible=false, ticklabelsvisible=false)
# colsize!(fig.layout, 1, Aspect(1, 1.0))

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

alert("'T is ier allemolle gedoan")


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
    plot_name_1 = "osc_α=$(αs[1]).pdf"
    save(plot_name_1, fig)
    plot_name_2 = "coex_α=$(αs[1]).pdf"
    save(plot_name_2, fig2)
    plot_name_3 = "coex_osc_α=$(αs[1]).pdf"
    save(plot_name_3, fig3)
    println("saved $plot_name_3")
end

fig3

# if isdir(datadir("b2basics/type3")) == false
#     mkpath(datadir("b2basics/type3"))
# end

# jldsave(datadir("b2basics/type3/oscillation_matrix_presentatie.jld2"); convergence_matrix, coexistence_matrix, ω_range, β_range)

# dfje = collect_results!(datadir("b2basics/type3"))




# save_dir = "C:/Users/tvgiel/OneDrive - UGent/Documents/LaTeX/Thomas VAN GIEL, Higher order interactions/Natasha_model/Figures/back2basics/ABC_ABD/" 
# save_name = save_dir * "coexistence_matrix_ωβ2.svg"

# save(save_name, fig)