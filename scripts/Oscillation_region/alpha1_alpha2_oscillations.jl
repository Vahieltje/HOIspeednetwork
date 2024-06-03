using DrWatson
@quickactivate "Nataha community model"
include(srcdir("Oscillation_region/bifurcation_alpha1_alpha2_functions.jl"))

using Distributed
using Alert
CairoMakie.activate!()


##############
# parameters #
##############

# the range of ω's and β's to consider per pixel of the final image. These whole ranges are considered for each pixel
ω_range = -3:0.4:2
β_range = -80:4:0


# The range of α values. One combination of α12 and α23 is one pixel
α12_range = 0:0.03:3
α23_range = 0:0.03:3

HOI_types = (asymmetric_HOI_ecosystem, symmetric_HOI_ecosystem)
ijk_species_list = [(1, 2, 3), (1, 3, 2), (2, 3, 1)]
intransitive = false
boosted_αs = ["α12", "α13", "α23"]

ijk_species = ijk_species_list[2]         # (i,j,k) species. k modifies the effect of j on i 
HOI_type = asymmetric_HOI_ecosystem       # Should be either asymmetric_HOI_ecosystem or symmetric_HOI_ecosystem, depending on if the modification is symmetric or asymmetric
boosted_α = boosted_αs[2]   # which α to be non-identical to the other two.
                            # α12 = αAB, α13 = αAC, α23 = αBC

savefig = true



α12_α23_file = filename_α(HOI_type, ijk_species, boosted_α, intransitive)
println("α12_α23_file = $α12_α23_file")

println("α12_α23_file = $α12_α23_file")

# check if the data is already calculated and saved
# If the data doesn't exist yet, calculate it and save it
if !isfile(α12_α23_file)
  α12_α23_oscillations = zeros(length(α12_range), length(α23_range))

  # from https://github.com/Arpeggeo/julia-distributed-computing. This loop is used for parallel computing
  status = pmap(1:length(α12_range)) do i

      α12 = α12_range[i]
      for (j, α23) in enumerate(α23_range)
        # set the correct α values, depending on the non-identical α choice
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

        # Check if the data of one pixel is already calculated and saved. 
        oscillation_file = filepath_αs(; data_type="oscillation", αs, HOI_type, ijk_species, intransitive)
      coexistence_file = filepath_αs(; data_type="coexistence", αs, HOI_type, ijk_species, intransitive)
        # if the data is not saved, calculate it and save it
        if isfile(oscillation_file) == false
          # calculate the data
          convergence_matrix, coexistence_matrix = make_bifurcation_matrix(ω_range, β_range, ijk_species; αs, HOI_type, intransitive)
          # save the data
          CSV.write(oscillation_file, Tables.table(convergence_matrix), writeheader=false)
          CSV.write(coexistence_file, Tables.table(coexistence_matrix), writeheader=false)
        else
          convergence_matrix = Matrix(CSV.read(oscillation_file, DataFrame))
        end

        # change the value of the pixel in the matrix to the fraction of systems that oscillate
        α12_α23_oscillations[i, j] = (length(convergence_matrix) - sum(convergence_matrix)) / length(convergence_matrix)
        println(" $(α12_α23_oscillations[i, j])")
      end
  end

  # write the full matrix to a file
  print("writing to: ")
  println(α12_α23_file)
  CSV.write(α12_α23_file, Tables.table(α12_α23_oscillations), writeheader=false)

# If the data is already calculated and saved, load it
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