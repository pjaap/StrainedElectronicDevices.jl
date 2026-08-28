module Marcogliese_et_al

using StrainedElectronicDevices

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
const cell_region_VS = 1
const cell_region_elec = 2
const cell_region_stressor = 3
const cell_region_homogen = 4

# boundary region assignments
const boundary_region_left = 5
const boundary_region_right = 3
const boundary_region_front = 2
const boundary_region_back = 4
const boundary_region_top = 6
const boundary_region_bottom = 1
const boundary_region_default = 7

function make_grid(;
        w_b = 200_000.0, # width plate
        θ = 54.7, # angle diag wall in degrees
        d_x = 525_000.0, # height of etched hole
        t_m = 1400.0, # silicon layer thickness
        t_e = 50.0, # electronic layer thickness
        t_s = 500.0, # height of stressor
        w_s = 5000.0, # width of stressor
        l_s = 10_000.0, # length of stressor
        d_s = 150.0, # distance between stressors
        w_r = 300_000.0, # dist to boundary
        w_i = 20_000.0, # width/length of region of interest
        kwargs...
    )

    @assert w_i ≥ max(2w_s + d_s, l_s) "stressor must fit into the region of interest"

    δ = d_x / tan(θ / 360 * 2π) |> abs
    w_m = w_b + 2δ

    # heights (z axis)
    h1 = 0.0 #-(d_x + t_m + t_e) # then: stressors are at level z = 0
    h3 = h1 + d_x
    h4 = h3 + t_m
    h5 = h4 + t_e
    h6 = h5 + t_s

    builder = SimplexGridBuilder(; Generator = TetGen)

    # outer bottom
    r1 = w_m / 2
    r2 = r1 + w_r
    p00 = point!(builder, -r2, -r2, h1)
    p10 = point!(builder, r2, -r2, h1)
    p01 = point!(builder, -r2, r2, h1)
    p11 = point!(builder, r2, r2, h1)

    # inner bottom
    q00 = point!(builder, -r1, -r1, h1)
    q10 = point!(builder, r1, -r1, h1)
    q01 = point!(builder, -r1, r1, h1)
    q11 = point!(builder, r1, r1, h1)

    # bottom corners of plate
    r3 = w_b / 2
    r00 = point!(builder, -r3, -r3, h3)
    r10 = point!(builder, r3, -r3, h3)
    r01 = point!(builder, -r3, r3, h3)
    r11 = point!(builder, r3, r3, h3)

    # region of interest
    r4 = w_i / 2
    i00 = point!(builder, -r4, -r4, h3)
    i10 = point!(builder, r4, -r4, h3)
    i01 = point!(builder, -r4, r4, h3)
    i11 = point!(builder, r4, r4, h3)

    j00 = point!(builder, -r4, -r4, h4)
    j10 = point!(builder, r4, -r4, h4)
    j01 = point!(builder, -r4, r4, h4)
    j11 = point!(builder, r4, r4, h4)

    k00 = point!(builder, -r4, -r4, h5)
    k10 = point!(builder, r4, -r4, h5)
    k01 = point!(builder, -r4, r4, h5)
    k11 = point!(builder, r4, r4, h5)

    # top corner points
    t00 = point!(builder, -r2, -r2, h5)
    t10 = point!(builder, r2, -r2, h5)
    t01 = point!(builder, -r2, r2, h5)
    t11 = point!(builder, r2, r2, h5)

    # stressor corners
    r5 = d_s / 2
    r6 = r5 + w_s
    r7 = l_s / 2
    s00 = point!(builder, -r6, -r7, h5)
    s10 = point!(builder, -r5, -r7, h5)
    s20 = point!(builder, r5, -r7, h5)
    s30 = point!(builder, r6, -r7, h5)

    s01 = point!(builder, -r6, r7, h5)
    s11 = point!(builder, -r5, r7, h5)
    s21 = point!(builder, r5, r7, h5)
    s31 = point!(builder, r6, r7, h5)

    # upper stressor corners
    u00 = point!(builder, -r6, -r7, h6)
    u10 = point!(builder, -r5, -r7, h6)
    u20 = point!(builder, r5, -r7, h6)
    u30 = point!(builder, r6, -r7, h6)

    u01 = point!(builder, -r6, r7, h6)
    u11 = point!(builder, -r5, r7, h6)
    u21 = point!(builder, r5, r7, h6)
    u31 = point!(builder, r6, r7, h6)

    # bottom
    facetregion!(builder, boundary_region_bottom)
    facet!(builder, [p00, q00, q01, q11, p11, p01])
    facet!(builder, [p00, p10, p11, q11, q10, q00])

    # bottom ramp
    facetregion!(builder, boundary_region_default)
    facet!(builder, q00, r00, r01, q01)
    facet!(builder, q00, q10, r10, r00)
    facet!(builder, q10, q11, r11, r10)
    facet!(builder, q01, q11, r11, r01)

    # bottom of plate
    facet!(builder, i00, i10, i11, i01)
    facet!(builder, i00, r00, r01, i01)
    facet!(builder, i00, i10, r10, r00)
    facet!(builder, i10, i11, r11, r10)
    facet!(builder, r01, r11, i11, i01)

    # side facets
    facetregion!(builder, boundary_region_left)
    facet!(builder, p00, p01, t01, t00)

    facetregion!(builder, boundary_region_right)
    facet!(builder, p10, p11, t11, t10)

    facetregion!(builder, boundary_region_front)
    facet!(builder, p00, p10, t10, t00)

    facetregion!(builder, boundary_region_back)
    facet!(builder, p01, p11, t11, t01)

    # top facets
    facetregion!(builder, boundary_region_top)
    facet!(builder, [t00, k00, k10, k11, t11, t10])
    facet!(builder, [t00, k00, k01, k11, t11, t01])
    facet!(builder, [k00, k10, k11, s31, s30, s20, s10, s00])
    facet!(builder, [k00, s00, s01, s11, s21, s31, k11, k01])
    facet!(builder, s10, s20, s21, s11)


    # internal facets
    facetregion!(builder, 1)

    facet!(builder, j00, j01, j11, j10)

    facet!(builder, i00, i10, j10, j00)
    facet!(builder, i10, i11, j11, j10)
    facet!(builder, i01, i11, j11, j01)
    facet!(builder, i00, i01, j01, j00)

    facetregion!(builder, 2)

    facet!(builder, j00, j10, k10, k00)
    facet!(builder, j10, j11, k11, k10)
    facet!(builder, j01, j11, k11, k01)
    facet!(builder, j00, j01, k01, k00)

    # # stressors
    facetregion!(builder, 8)
    facet!(builder, s00, s10, s11, s01)
    facet!(builder, s20, s30, s31, s21)

    facet!(builder, s00, s10, u10, u00)
    facet!(builder, s10, s11, u11, u10)
    facet!(builder, s01, s11, u11, u01)
    facet!(builder, s00, s01, u01, u00)
    facet!(builder, u00, u10, u11, u01)

    facet!(builder, s20, s30, u30, u20)
    facet!(builder, s30, s31, u31, u30)
    facet!(builder, s21, s31, u31, u21)
    facet!(builder, s20, s21, u21, u20)
    facet!(builder, u20, u30, u31, u21)

    # options!(builder; radius_edge_ratio = 2)

    cellregion!(builder, cell_region_homogen)
    maxvolume!(builder, Inf64)
    regionpoint!(builder, -(r1 + r2) / 2, 0, (h1 + h3) / 2)

    cellregion!(builder, cell_region_VS)
    maxvolume!(builder, Inf64)
    regionpoint!(builder, 0, 0, (h3 + h4) / 2)

    cellregion!(builder, cell_region_elec)
    maxvolume!(builder, Inf64)
    regionpoint!(builder, 0, 0, (h4 + h5) / 2)

    cellregion!(builder, cell_region_stressor)
    maxvolume!(builder, Inf64)
    regionpoint!(builder, r5 + w_s / 2, 0, h5 + t_s / 2)
    regionpoint!(builder, -r5 - w_s / 2, 0, h5 + t_s / 2)


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
        stress = 560.0,
        kwargs...
    )

    xgrid = uniform_refine(make_grid(kwargs...), nref)

    # workaround https://github.com/WIAS-PDELib/ExtendableGrids.jl/issues/136
    z_shift = -526450.0
    xgrid[Coordinates][3, :] .+= z_shift

    npart = 9 * Threads.nthreads()
    xgrid = partition(xgrid, PlainMetisPartitioning(; npart))
    @info "done partitioning the grid into $npart parts with partitions per color = $(num_partitions_per_color(xgrid))"

    materials = material_vector(4)
    materials[cell_region_VS] = SiGe(0.3)
    materials[cell_region_elec] = Al()
    materials[cell_region_stressor] = TiN(:A)
    materials[cell_region_homogen] = SiGe(0.5)


    # unit matrix in Voigt notation
    Iᵥ = @SArray [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

    thermal_strain = [
        cell_region_stressor => Iᵥ * stress,
    ]

    # create the electronic device
    device = Device(xgrid, materials; thermal_strain)

    # create the linear elasticity problem
    elasticity_problem = create_linear_elasticity_problem(
        device;
        dirichlet_boundary = [boundary_region_bottom => 0.0]
    )

    sol_elasticity = simulate_elasticity(elasticity_problem, xgrid; order = order_displacement)

    return sol_elasticity, device
end


function plot(
        sol_elasticity,
        device;
        kwargs...
    )

    ## displacement is the first component of the solution
    displacement = sol_elasticity.tags[1]

    # the grid with all adjacencies removed (for discontinuous plotting)
    xgrid = explode(sol_elasticity[displacement].FES.xgrid)

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
    lazy_interpolate!(strain_func[1], sol_elasticity, [εV(displacement, 1.0)], postprocess = add_pre_strain_kernel!, use_cellparents = true)

    # create a strain FE function
    FES_displacement = FESpace{H1P1(3)}(xgrid)
    displacement_func = FEVector(FES_displacement)
    lazy_interpolate!(displacement_func[1], sol_elasticity, use_cellparents = true)

    vis = GridVisualizer(Plotter = GLMakie, size = (1500, 1200), layout = (2, 3), show = false)
    strain_vals = nodevalues(strain_func[1])
    @views scalarplot!(vis[1, 1], xgrid, strain_vals[1, :], title = "ε₁₁", slice = :z => -10.0)
    @views scalarplot!(vis[1, 2], xgrid, strain_vals[2, :], title = "ε₂₂", slice = :z => -10.0)
    @views scalarplot!(vis[1, 3], xgrid, strain_vals[3, :], title = "ε₃₃", slice = :z => -10.0)
    @views scalarplot!(vis[2, 1], xgrid, strain_vals[4, :], title = "ε₂₃", slice = :z => -10.0)
    @views scalarplot!(vis[2, 2], xgrid, strain_vals[5, :], title = "ε₁₃", slice = :z => -10.0)
    @views scalarplot!(vis[2, 3], xgrid, strain_vals[6, :], title = "ε₁₂", slice = :z => -10.0)
    reveal(vis)

    # export the nodevalues to VTK
    writeVTK(
        "Marcogliese_et_al.vtu",
        xgrid;
        compress = true,
        :cell_regions => xgrid[CellRegions],
        :displacement => nodevalues(displacement_func[1]),
        :strain => nodevalues(strain_func[1]),
    )

    return nothing

end

function main(; kwargs...)
    result = simulate(; kwargs...)
    plot(result...; kwargs...)
    return result
end

end # module
