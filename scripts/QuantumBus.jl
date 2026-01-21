module QuantumBus

using StrainedElectronicDevices

using StaticArrays
using ExtendableFEM
using ExtendableGrids
using GridVisualize
using JLD2: load, save
# using PythonCall

# the default linear solver
using Pardiso
using LinearSolve: PardisoJL

const dim = 3

# grid region assignments
const cell_region_Si = 1
const cell_region_SiGe = 2
const cell_region_SiO₂ = 3
const cell_region_TiN = 4

# boundary region assignments
const boundary_region_left = 2
const boundary_region_right = 3
const boundary_region_bottom = 4

function process_grid(xgrid)

    # post-process the boundary regions: define "left", "right" and "bottom" region
    grid_x_range = @views extrema(xgrid[Coordinates][1, :])
    grid_y_range = @views extrema(xgrid[Coordinates][2, :])
    grid_z_range = @views extrema(xgrid[Coordinates][3, :])

    # reset all face regions
    bfacemask!(
        xgrid,
        (grid_x_range[1], grid_y_range[1], grid_z_range[1]),
        (grid_x_range[2], grid_y_range[2], grid_z_range[2]),
        1
    )

    # set left and right periodic boundary regions
    bfacemask!(
        xgrid,
        (grid_x_range[1], grid_y_range[1], grid_z_range[1]),
        (grid_x_range[1], grid_y_range[2], grid_z_range[2]),
        boundary_region_left
    )

    bfacemask!(
        xgrid,
        (grid_x_range[2], grid_y_range[1], grid_z_range[1]),
        (grid_x_range[2], grid_y_range[2], grid_z_range[2]),
        boundary_region_right
    )

    # set bottom region
    bfacemask!(
        xgrid,
        (grid_x_range[1], grid_y_range[1], grid_z_range[1]),
        (grid_x_range[2], grid_y_range[2], grid_z_range[1]),
        boundary_region_bottom
    )

    return xgrid
end


function simulate(;
        linear_solver = PardisoJL,
        TiN_mode = :A, # choose :A or :B
        grid_variant = :coarse, # choose :coarse or :fine
        σ_0 = -2.6,
        periodic = true,
        kwargs...
    )

    materials = material_vector(4)
    materials[cell_region_Si] = Si()
    materials[cell_region_SiGe] = SiGe(0.34)
    materials[cell_region_SiO₂] = SiO₂()
    materials[cell_region_TiN] = TiN(TiN_mode)

    # unit matrix in Voigt notation
    Iᵥ = @SArray [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

    # in-plane (x-y) stress unit matrix in Voigt notation
    Jᵥ = @SArray [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

    # external pre-stress at the TiN claviers
    pre_stress = [
        cell_region_TiN => σ_0 * Iᵥ,
    ]

    ## read the grid from a file
    xgrid = load(StrainedElectronicDevices.gridsdir("qubus_one_contact_$(grid_variant).jld2"))["grid"]
    @info "grid loaded"

    # process the grid (add additional boundary regions)
    xgrid = process_grid(xgrid)

    # lift the grid coordinates to nm (otherwise the cell volumes become to numerically zero)
    xgrid[Coordinates] *= 1.0e6

    # create the electronic device
    device = Device(xgrid, materials; pre_stress)

    # create the linear elasticity problem
    elasticity_problem = create_linear_elasticity_problem(
        device;
        dirichlet_boundary = [boundary_region_bottom => 0.0],
        periodic_coupling = periodic ? [boundary_region_left => boundary_region_right] : []
    )

    # second order finite element space
    FES = FESpace{H1P2{3, 3}}(xgrid)

    #solve
    sol = ExtendableFEM.solve(
        elasticity_problem,
        FES;
        method_linear = linear_solver
    )

    return sol, device
end


function plot(
        sol,
        device;
        Plotter = nothing, # include a plotter in your global environment: GLMakie, PythonPlot...
        kwargs...
    )

    ## displacement is the first component of the solution
    displacement = sol.tags[1]
    xgrid = sol[displacement].FES.xgrid
    @showtime xgrid = explode(xgrid)


    # extract pre-strains (from pre-stress)
    pre_strain = [ material_tensor \ pre_stress for (material_tensor, pre_stress) in zip(device.material_tensors, device.pre_stress) ]

    # create a strain FE function
    FES_strain = FESpace{H1P1(6)}(xgrid)

    # post process interpolator
    function add_pre_strain_kernel!(result, input, qpinfo)
        @. result = input + pre_strain[qpinfo.region]
        return nothing
    end

    strain_func = FEVector(FES_strain)
    lazy_interpolate!(strain_func[1], sol, [εV(displacement, 1.0)], postprocess = add_pre_strain_kernel!, use_cellparents = true)

    # create a strain FE function
    FES_displacement = FESpace{H1P1(3)}(xgrid)
    displacement_func = FEVector(FES_displacement)
    lazy_interpolate!(displacement_func[1], sol, use_cellparents = true)

    strain_values = nodevalues(strain_func[1])
    displacement_values = nodevalues(displacement_func[1])

    writeVTK(
        "QuantumBusResult.vtu",
        xgrid;
        :displacement => displacement_values,
        :strain => strain_values,
    )

    gridplot(xgrid; Plotter, yplanes = [0.00000039], scene3d = :LScene)


    return nothing

end

function main(; kwargs...)
    result = simulate(; kwargs...)
    plot(result...; kwargs...)
    return result
end

end # module
