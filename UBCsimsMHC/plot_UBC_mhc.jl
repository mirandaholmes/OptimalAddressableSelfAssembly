#Created Summer 2024 by Sushrut Tadwalkar. Makes an animation of a simulation of "UBC" assembly with diffusion. This particular version plots specific frames with a light grey background and the monomers having no borders.

using Plots #For the plotting
using ColorSchemes 
using Colors #For the UBC colour scheme.
using DelimitedFiles #To read data files.
#using Measures

#Reads the trajectory data file, and compiles it into a matrix of integers.
ifanimate = false
#num=4
#framemax = 129
#fpl = [framemax]
 num = 1
 framemax = 129   # max number of frames. set to "0" to plot all frames in simulation. 
 fpl = []#[1,2,3,5,7,17,25,97,98]    # which frames to plot and save (empty for none)
#num = 4
#framemax = 62
#fpl = [framemax]
filename = "Data/ubc_mhc_"*string(num)*"_traj.txt"
filepath = joinpath(@__DIR__, filename) #Creates a relative path.
system_sizes = readdlm(filepath)[2,:] #Reads the system sizes.
coordinates = readdlm(filepath, skipstart = 3)[:,[2,3]] #Data in the form of a matrix of integers.

#System sizes.
nParticles = system_sizes[1] #Number of monomers.
boxLength = system_sizes[2] #Length of the box.
n0 = system_sizes[3] #Number of monomers in the base cell.
l0 = sqrt(n0) #Side length of each monomer (will be defined differently for larger systems.)

#Making the UBC color map.
UBC_colors = [RGB(12/255,35/255,68/255),RGB(0,85/255,183/255),RGB(0,119/255,200/255),RGB(0,167/255,225/255),RGB(64/255,180/255,229/255),RGB(110/255,196/255,232/255),RGB(1,1,1)] #Making the color array.

UBC_colorscheme = ColorScheme(UBC_colors) #Turning the array into a color scheme.

UBC_colormap = get.(Ref(UBC_colorscheme), range(0, 1, length=256)) #Making a color map from the color scheme.
#------------------------------------------------------------------------------#

#Following code sets up important functions that are used (directly/indirectly) in plotting.

N = size(coordinates)[1]

square_vertices = Vector{Shape}()

#Drawing squares and storing them as vertices.
function make_vertices_array(args...)
    for j in 1:N
        if typeof(coordinates[j,1]) == Float64
            push!(square_vertices, Shape(coordinates[j,1]  .+ [0,1,1,0], coordinates[j,2] .+ [0,0,1,1]))
            # this one is for num=1, fpl=97
            #push!(square_vertices, Shape(mod(coordinates[j,1] + 80.0, boxLength) .+ [0,1,1,0], coordinates[j,2] .+ [0,0,1,1]))
            # this one is for num=1, fpl=98
            #push!(square_vertices, Shape(mod(coordinates[j,1] + 46.0, boxLength) .+ [0,1,1,0], coordinates[j,2] .+ [0,0,1,1]))
            end
    end
end

make_vertices_array()

n = Int.(size(square_vertices)[1]/nParticles)

squares = Vector{Vector{Shape}}() #Array to store the frames that will be plotted.

#Separating the squares by frames.
function make_plot_array(args...)
    for i in 0:n-1
        squares0 = Vector{Shape}()
        for j in 1:nParticles
            push!(squares0, square_vertices[j+nParticles*i])
        end
        push!(squares, squares0)
    end
end

make_plot_array()

#Colour assigning function.
function assign_color(index = 1)
    residue = mod(index, n0)
    if residue === 0
        residue = n0
    end
    color_code = Int.(ceil(residue*256/n0))
    UBC_colormap[color_code]
end

colours0 = Vector{RGB}()

function make_color_vector(args...)
    for i in 1:nParticles
        push!(colours0, assign_color(i))
    end
end

make_color_vector()

colours = permutedims(colours0) #Plotting requires a row vector.

#------------------------------------------------------------------------------#

if framemax != 0
    frames = size(squares)[1]   # number of frames
end

# Making the animation.
if ifanimate
    plot_animation= @animate for i in 1:frames
    plot(squares[i], xlim = (0,boxLength), ylim = (0,boxLength), title=string(i-1)*" x 2e4 MC sweeps", aspect_ratio=:equal, legend = false, fillcolor = colours, axis=false, framestyle=:none,
          linewidth = 0, linecolor =  colours, background_color =:lightgrey,size=(500, 500))
    end
    filename2 = "Plots/mv_UBC_"*string(num)*".mp4"
    filepath2 = joinpath(@__DIR__, filename2) #Creates a relative path.
    gif(plot_animation, filepath2, fps=5)
end

#Plots a specific frame; generally used for plotting last frame, for which the argument passed in should be "frames".
function psf(fn = 1)
    plot(squares[fn], xlim = (0,boxLength), ylim = (0,boxLength), title=string(fn-1)*" x 2e4 MC sweeps", aspect_ratio=:equal, legend = false, fillcolor = colours, axis=false, framestyle=:none, 
    linewidth = 0, linecolor =  colours, background_color =:lightgrey,size=(500, 500))
end





for j=1:length(fpl)
    jf = fpl[j]
    psf(fpl[j])
    filepath = joinpath(@__DIR__, "Plots/pl_UBC_"*string(num)*"_t"*string(jf-1)*".pdf") #Creates a relative path.
    savefig(filepath)
end

# psf(8) #This one clips out because of boundary conditions, so replace "coordinates[j,1] .+ [0,1,1,0]" with "mod(coordinates[j,1] + 50.0, boxLength) .+ [0,1,1,0]" in make_vertices_array
# filename10 = "Plots/nb_UBC_8.pdf"
# filepath10 = joinpath(@__DIR__, filename10) #Creates a relative path.
# savefig(filepath10)
