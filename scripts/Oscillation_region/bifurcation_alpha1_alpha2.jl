using DrWatson
@quickactivate "Nataha community model"
include(srcdir("Oscillation_region/bifurcation_alpha1_alpha2_functions.jl"))

using Distributed
using Alert
CairoMakie.activate!()
##############
# parameters #
##############


ω_range = -3:0.1:2
β_range = -80:1:0
# β_range = 0:0.1:40
α12_range = 0:0.03:3
α23_range = 0:0.03:3


HOI_types = (single_HOI_ecosystem, symmetric_HOI_ecosystem)
ijk_species_list = [(1, 2, 3), (1, 3, 2), (2, 3, 1)]
intransitive = false
superiority = (1, 2, 3)
boosted_αs = ["α12", "α13", "α23"]

ijk_species = ijk_species_list[2]
HOI_type = single_HOI_ecosystem
boosted_α = boosted_αs[2]




profiling = false
savefig = true

α12_α23_file = filename_α(HOI_type, ijk_species, boosted_α, intransitive)
println("α12_α23_file = $α12_α23_file")

println("α12_α23_file = $α12_α23_file")

if !isfile(α12_α23_file) || profiling
  α12_α23_oscillations = zeros(length(α12_range), length(α23_range))

  # from https://github.com/Arpeggeo/julia-distributed-computing
  status = pmap(1:length(α12_range)) do i
    #try
      α12 = α12_range[i]
      for (j, α23) in enumerate(α23_range)
        # vergeet naam α12_α23_file ook niet aan te passen
        if boosted_α == "α12"
          αs = [α12, α23, α23]
        elseif boosted_α == "α13"
          αs = [α23, α12, α23]
        elseif boosted_α == "α23"
          αs = [α23, α23, α12]
        else
          error("boosted_α should be either 'α12', 'α13' or 'α23'")
        end
        

        print("αs = $αs, HOI_type = $HOI_type, ijk_species = $ijk_species :")
        oscillation_file = filepath_omega_beta(; data_type="oscillation", αs, HOI_type, ijk_species, intransitive, superiority)
        coexistence_file = filepath_omega_beta(; data_type="coexistence", αs, HOI_type, ijk_species, intransitive, superiority)

        if isfile(oscillation_file) == false || profiling
          convergence_matrix, coexistence_matrix = make_bifurcation_matrix(ω_range, β_range, ijk_species; αs, HOI_type, intransitive, superiority)
          CSV.write(oscillation_file, Tables.table(convergence_matrix), writeheader=false)
          CSV.write(coexistence_file, Tables.table(coexistence_matrix), writeheader=false)
        else
          convergence_matrix = Matrix(CSV.read(oscillation_file, DataFrame))
        end

        α12_α23_oscillations[i, j] = (length(convergence_matrix) - sum(convergence_matrix)) / length(convergence_matrix)
        println(" $(α12_α23_oscillations[i, j])")
      end
    # catch e
    #   false
    # end
  end
  print("writing to: ")
  println(α12_α23_file)
  CSV.write(α12_α23_file, Tables.table(α12_α23_oscillations), writeheader=false)
else
  α12_α23_oscillations = Matrix(CSV.read(α12_α23_file, DataFrame))
end

# get_new_data = true
# save_new_data = false & get_new_data

############
# Plotting #
############

colormap_oscillation = [ugent_blue, ugent_yellow]
colormap_oscillation_inverted = [colormap_oscillation[2], colormap_oscillation[1]]
# Oscillation heatmap
println("plotting heatmap")
fig = Figure()
ax = Axis(fig[1, 1], xlabel="α₁₂", ylabel="α₁₃ and α₂₃", aspect=1, #title="HOI_type = $HOI_type, ijk_species = $ijk_species"
)
cmap = cgrad(:viridis, categorical=false)
hm = Makie.heatmap!(ax, α12_range, α23_range, α12_α23_oscillations, colormap=cmap, colorrange=(0, 0.15))


subgrid = GridLayout(fig[1, 2], tellheight=false)
# Label(subgrid[1, 1], "Some combinations\noscillating")
# Label(subgrid[3, 1], "No combinations\noscillating")
cb = Colorbar(
  subgrid[1, 1],
  hm
  # colormap=cgrad(colormap_oscillation, 2, categorical=false),
  # ticksvisible=false, ticklabelsvisible=false
)


colsize!(fig.layout, 1, Aspect(1, 1.0))

resize_to_layout!(fig)
fig

if savefig
  save_folder = "paper/"
  plot_name = "temp_α_α$(ijk_species).pdf"
  save(plotsdir(plot_name), fig)
end

alert("'T is ier allemolle gedoan")

fig