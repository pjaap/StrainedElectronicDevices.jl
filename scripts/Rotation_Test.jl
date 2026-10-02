module Rotation_Test

using StrainedElectronicDevices

using LinearAlgebra
using StaticArrays
using ExtendableFEM
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using JLD2: load, save
using LinearAlgebra: I, Diagonal, lu, diag, Symmetric
using SparseArrays: sparse
using TetGen
using SimplexGridFactory
using Metis

using UnicodePlots
using GLMakie

# the default linear solver
using Pardiso
using Krylov
using LinearSolve

const dim = 3

# grid region assignments
const cell_region_bulk = 1
const cell_region_stressor1 = 2
const cell_region_stressor2 = 3

# boundary region assignments
const boundary_region_left = 5
const boundary_region_right = 3
const boundary_region_front = 2
const boundary_region_back = 4
const boundary_region_top = 6
const boundary_region_bottom = 1
const boundary_region_default = 7

function make_grid(;
        w_b = 20_000.0, # width bulk plate
        h_b = 1_000.0, # height of bulk
        h_s = 500.0, # height of stressor
        w_s = 5_000.0, # width of stressor
        d_s = 150.0, # distance between stressors
        max_vol = 1.0e5, # max grid element volume
        fine_vol = max_vol / 1.0e4, # fine grained volume
        slice_thickness = 100.0,
        kwargs...
    )

    # @assert w_b ≥ max(2w_s + d_s, l_s) "stressor must fit into the bulk plate"


    # heights (z axis)
    h1 = -h_b
    h2 = 0.0
    h3 = h_s

    builder = SimplexGridBuilder(; Generator = TetGen)

    # outer bottom
    r2 = w_b / 2
    r3 = slice_thickness / 2
    p00 = point!(builder, -r2, -r3, h1)
    p10 = point!(builder, r2, -r3, h1)
    p01 = point!(builder, -r2, r3, h1)
    p11 = point!(builder, r2, r3, h1)


    # top corner points
    t00 = point!(builder, -r2, -r3, h2)
    t10 = point!(builder, r2, -r3, h2)
    t01 = point!(builder, -r2, r3, h2)
    t11 = point!(builder, r2, r3, h2)

    # stressor corners
    r5 = d_s / 2
    r6 = r5 + w_s
    s00 = point!(builder, -r6, -r3, h2)
    s10 = point!(builder, -r5, -r3, h2)
    s20 = point!(builder, r5, -r3, h2)
    s30 = point!(builder, r6, -r3, h2)

    s01 = point!(builder, -r6, r3, h2)
    s11 = point!(builder, -r5, r3, h2)
    s21 = point!(builder, r5, r3, h2)
    s31 = point!(builder, r6, r3, h2)

    # upper stressor corners
    u00 = point!(builder, -r6, -r3, h3)
    u10 = point!(builder, -r5, -r3, h3)
    u20 = point!(builder, r5, -r3, h3)
    u30 = point!(builder, r6, -r3, h3)

    u01 = point!(builder, -r6, r3, h3)
    u11 = point!(builder, -r5, r3, h3)
    u21 = point!(builder, r5, r3, h3)
    u31 = point!(builder, r6, r3, h3)

    # bottom
    facetregion!(builder, boundary_region_bottom)
    facet!(builder, p00, p01, p11, p10)

    # side facets
    facetregion!(builder, boundary_region_left)
    facet!(builder, p00, p01, t01, t00)

    facetregion!(builder, boundary_region_right)
    facet!(builder, p10, p11, t11, t10)

    facetregion!(builder, boundary_region_front)
    facet!(builder, [p00, p10, t10, s30, s20, s10, s00, t00])
    facet!(builder, s00, s10, u10, u00)
    facet!(builder, s20, s30, u30, u20)

    facetregion!(builder, boundary_region_back)
    facet!(builder, [p01, p11, t11, s31, s21, s11, s01, t01])
    facet!(builder, s01, s11, u11, u01)
    facet!(builder, s21, s31, u31, u21)

    # top facets
    facetregion!(builder, boundary_region_top)
    facet!(builder, t00, s00, s01, t01)
    facet!(builder, s30, t10, t11, s31)
    facet!(builder, s10, s20, s21, s11)

    facetregion!(builder, boundary_region_default)

    # stressors
    facet!(builder, s00, s10, s11, s01)
    facet!(builder, s20, s30, s31, s21)

    facet!(builder, s10, s11, u11, u10)
    facet!(builder, s00, s01, u01, u00)
    facet!(builder, u00, u10, u11, u01)

    facet!(builder, s30, s31, u31, u30)
    facet!(builder, s20, s21, u21, u20)
    facet!(builder, u20, u30, u31, u21)

    cellregion!(builder, cell_region_bulk)
    maxvolume!(builder, max_vol)
    regionpoint!(builder, 0, 0, (h1 + h2) / 2)

    cellregion!(builder, cell_region_stressor1)

    x_z_center = [0.0, -40.0]
    function unsuitable(p1, p2, p3, p4)
        vol = abs(det([p1 - p2 p1 - p3 p1 - p4])) / 2
        center = (p1 + p2 + p3 + p4) / 4
        dist = norm(center[[1, 3]] - x_z_center)

        desired_vol = dist > d_s ? max_vol : fine_vol
        return vol > desired_vol
    end
    options!(builder; unsuitable = unsuitable)

    regionpoint!(builder, -r5 - w_s / 2, 0, (h2 + h3) / 2)

    cellregion!(builder, cell_region_stressor2)
    maxvolume!(builder, max_vol)
    regionpoint!(builder, r5 + w_s / 2, 0, (h2 + h3) / 2)


    return simplexgrid(builder)
end


function simulate_elasticity(elasticity_problem_problem, xgrid; order)

    if order == 1
        FES = FESpace{H1P1{3}}(xgrid)
    elseif order == 2
        FES = FESpace{H1P2{3, 3}}(xgrid)
    else
        error("supported FE orders are 1 and 2.")
    end


    sol = ExtendableFEM.solve(
        elasticity_problem_problem,
        FES;
        parallel = true,
        verbosity = 2,
        method_linear = PardisoJL()
    )

    return sol
end

function simulate(;
        order_displacement = 1,
        nref = 0,
        stress_SiN = 1.0, # GPa
        T_final = 0.0, # K
        angle = 17π / 180,
        kwargs...
    )

    xgrid = uniform_refine(make_grid(; kwargs...), nref)

    npart = 9 * Threads.nthreads()
    xgrid = partition(xgrid, PlainMetisPartitioning(; npart))
    @info "done partitioning the grid into $npart parts with partitions per color = $(num_partitions_per_color(xgrid))"

    materials = material_vector(3)
    materials[cell_region_bulk] = Si()
    materials[cell_region_stressor1] = Si₃N₄()
    materials[cell_region_stressor2] = Si₃N₄()

    # x-yuUnit matrix in Voigt notation
    Jᵥ = @SArray [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

    pre_stress = [
        cell_region_stressor1 => Jᵥ * stress_SiN,
        cell_region_stressor2 => Jᵥ * stress_SiN,
    ]

    # rotate grid + materials
    rotation_matrix = z_rotation_matrix(angle)
    grid_rotated = deepcopy(xgrid)
    grid_rotated[Coordinates] = rotation_matrix * grid_rotated[Coordinates]

    device_default = Device(xgrid, materials; pre_stress)
    elasticity_problem_default = create_linear_elasticity_problem(
        device_default;
        dirichlet_boundary = [boundary_region_bottom => 0.0],
        periodic_coupling = [boundary_region_front => boundary_region_back]
    )
    sol_elasticity_default = simulate_elasticity(elasticity_problem_default, xgrid; order = order_displacement)


    device_rotated_grid = Device(grid_rotated, materials; pre_stress)
    elasticity_problem_rotated_grid = create_linear_elasticity_problem(
        device_rotated_grid;
        dirichlet_boundary = [boundary_region_bottom => 0.0],
        periodic_coupling = [boundary_region_front => boundary_region_back]
    )
    sol_elasticity_rotated_grid = simulate_elasticity(elasticity_problem_rotated_grid, grid_rotated; order = order_displacement)

    return sol_elasticity_default,
        sol_elasticity_rotated_grid,
        device_default,
        device_rotated_grid,
        rotation_matrix
end


function plot(
        sol_elasticity_default,
        sol_elasticity_rotated_grid,
        device_default,
        device_rotated_grid,
        rotation_matrix;
        kwargs...
    )

    xgrid = sol_elasticity_default[1].FES.xgrid
    grid_rotated = sol_elasticity_rotated_grid[1].FES.xgrid

    # the grid with all adjacencies removed (for discontinuous plotting)
    xgrid = explode(xgrid)
    grid_rotated = explode(grid_rotated)

    # extract pre-strains (from pre-stress)
    pre_strain1 = [ (pre_stress == zeros(6) ? pre_stress : material_tensor \ pre_stress) for (material_tensor, pre_stress) in zip(device_default.material_tensors, device_default.pre_stress) ]
    pre_strain3 = [ (pre_stress == zeros(6) ? pre_stress : material_tensor \ pre_stress) for (material_tensor, pre_stress) in zip(device_rotated_grid.material_tensors, device_rotated_grid.pre_stress) ]

    # create a strain FE function
    FES_strain = FESpace{H1P1(6)}(xgrid)
    FES_strain_rot = FESpace{H1P1(6)}(grid_rotated)

    # post process interpolator
    function make_prestrain_kernel(pre_strain)
        function add_pre_strain_kernel!(result, input, qpinfo)
            @. result = input + pre_strain[qpinfo.region]
            return nothing
        end
        return add_pre_strain_kernel!
    end

    strain_func1 = FEVector(FES_strain)
    strain_func3 = FEVector(FES_strain_rot)
    t1 = Threads.@spawn lazy_interpolate!(strain_func1[1], sol_elasticity_default, [εV(1, 1.0)], postprocess = make_prestrain_kernel(pre_strain1), use_cellparents = true)
    t2 = Threads.@spawn lazy_interpolate!(strain_func3[1], sol_elasticity_rotated_grid, [εV(1, 1.0)], postprocess = make_prestrain_kernel(pre_strain3), use_cellparents = true)
    fetch.([t1, t2])

    t1 = Threads.@spawn nodevalues(strain_func1[1])
    t2 = Threads.@spawn nodevalues(strain_func3[1])

    strain_vals1 = fetch(t1)
    strain_vals3 = fetch(t2)

    # # rotate back
    # s2v = StrainedElectronicDevices.strain2voigt
    # v2s = StrainedElectronicDevices.voigt2strain
    # R = rotation_matrix
    # for i in 1:size(strain_vals2, 2)
    #     strain_vals2[:, i] = s2v(R * v2s(strain_vals2[:, i]) * R')
    # end

    vis = GridVisualizer(Plotter = GLMakie, size = (1500, 1200), layout = (2, 3), show = false)
    @views scalarplot!(vis[1, 1], xgrid, strain_vals1[1, :], title = "def: ε₁₁", slice = :y => 0.0)
    @views scalarplot!(vis[1, 2], xgrid, strain_vals1[2, :], title = "def: ε₂₂", slice = :y => 0.0)
    @views scalarplot!(vis[1, 3], xgrid, strain_vals1[6, :], title = "def: ε₁₂", slice = :y => 0.0)

    @views scalarplot!(vis[2, 1], grid_rotated, strain_vals3[1, :], title = "rot. grid: ε₁₁", slice = :(x - y))
    @views scalarplot!(vis[2, 2], grid_rotated, strain_vals3[2, :], title = "rot. grid: ε₂₂", slice = :(x - y))
    @views scalarplot!(vis[2, 3], grid_rotated, strain_vals3[6, :], title = "rot. grid: ε₁₂", slice = :(x - y))

    reveal(vis)

    return nothing
end

function main(; kwargs...)
    result = simulate(; kwargs...)
    plot(result...; kwargs...)
    return result
end

end # module
