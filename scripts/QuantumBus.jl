module QuantumBus

using StrainedElectronicDevices

using StaticArrays
using ExtendableFEM
using ExtendableGrids
using GridVisualize
using JLD2: load, save
using LinearAlgebra: I, Diagonal
# using PythonCall

# the default linear solver
using Pardiso
using Krylov
using AMGCLWrap: AMGSolverAlgorithm, AMGPrecon
using LinearSolve: PardisoJL, KrylovJL_GMRES


using ILUZero: ILU0Precon

const dim = 3

# grid region assignments
const cell_region_Si = 1
const cell_region_SiGe = 2
const cell_region_SiO₂ = 3
const cell_region_TiN_side1 = 11 # side gates
const cell_region_TiN_side2 = 12
const cell_region_TiN_clav1 = 13 # clavier gates
const cell_region_TiN_clav2 = 14
const cell_region_TiN_clav3 = 15
const cell_region_TiN_clav4 = 16

# collect the TiN clavier regions
const cell_regions_TiN_clav = [
    cell_region_TiN_clav1,
    cell_region_TiN_clav2,
    cell_region_TiN_clav3,
    cell_region_TiN_clav4,
]

# collect the TiN side regions
const cell_regions_TiN_side = [
    cell_region_TiN_side1,
    cell_region_TiN_side2,
]

# boundary region assignments
const boundary_region_left = 11
const boundary_region_right = 12
const boundary_region_bottom = 13

function process_grid!(xgrid)

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

    return nothing
end


function simulate(;
        linear_solver = KrylovJL_GMRES(rtol = 1.0e-15, verbose = 10, precs = (A, p) -> (AMGPrecon(A), I)),
        TiN_mode = :A, # choose :A or :B
        grid_variant = :coarse, # choose :coarse or :fine
        σ_0 = -2.6,
        Plotter = nothing,
        periodic = true,
        order = 1,
    )

    materials = material_vector(16)
    materials[cell_region_Si] = Si()
    materials[cell_region_SiGe] = SiGe(0.34)
    materials[cell_region_SiO₂] = SiO₂()

    for cell_region in vcat(cell_regions_TiN_clav, cell_regions_TiN_side)
        materials[cell_region] = TiN(TiN_mode)
    end

    # in-plane (x-y) stress unit matrix in Voigt notation
    Jᵥ = @SArray [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

    # external pre-stress at the TiN claviers
    pre_stress = [
        cell_region => σ_0 * Jᵥ for cell_region in cell_regions_TiN_clav
    ]

    ## read the grid from a file
    xgrid = load(StrainedElectronicDevices.gridsdir("qubus_one_contact_$(grid_variant).jld2"))["grid"]
    @info "grid loaded"

    # process the grid (add additional boundary regions)
    process_grid!(xgrid)

    # lift the grid coordinates to nm (otherwise the cell volumes become to numerically zero)
    # and remove cached values from the grid (void losing integrity of the grid)
    xgrid[Coordinates] *= 1.0e9
    trim!(xgrid)

    # create the electronic device
    device = Device(xgrid, materials; pre_stress)

    # create the linear elasticity problem
    elasticity_problem = create_linear_elasticity_problem(
        device;
        dirichlet_boundary = [boundary_region_bottom => 0.0],
        periodic_coupling = periodic ? [boundary_region_left => boundary_region_right] : []
    )

    # second order finite element space
    if order == 1
        FES = FESpace{H1P1{3}}(xgrid) # only for local testing
    elseif order == 2
        FES = FESpace{H1P2{3, 3}}(xgrid)
    else
        error("supported FE orders are 1 and 2.")
    end

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

    # the grid with all adjacencies removed (for discontinuous plotting)
    xgrid = explode(sol[displacement].FES.xgrid)

    # extract pre-strains (from pre-stress)
    pre_strain = [ (pre_stress == zeros(6) ? pre_stress : material_tensor \ pre_stress) for (material_tensor, pre_stress) in zip(device.material_tensors, device.pre_stress) ]

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

    # export the nodevalues to VTK
    writeVTK(
        "QuantumBusResult.vtu",
        xgrid;
        compress = true,
        :displacement => nodevalues(displacement_func[1]),
        :strain => nodevalues(strain_func[1]),
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
