"""
    create_linear_elasticity_problem(device::Device; dirichlet_boundary = (), periodic_coupling::Union{Nothing, PeriodicCoupling} = nothing)

    This creates an ExtendableFEM compatible linear elasticity problem based on a given Device.

    Dirichlet data should be provided by a vector of pairs `[ region1 => value1, region2 => value2, ... ]``
    Periodic coupling should be provided by a vector of pairs `[] source_region1 => target_region1, ...]`

    The assembled PDE is reads

    div(σ + σ*) = 0              in Ω
              σ = Cε(u)          in Ω
           ε(u) = ½( ∇u + ∇uᵀ )  in Ω
              u = uᴰ             in Ωᴰ

    where
    - σ* is the total pre-stress (resulting from lattice mismatch, thermal strain and given external pre-stresses/strains)
    - uᴰ is given Dirichlet boundary data
    - C is the material tensor (location dependent)
"""
function create_linear_elasticity_problem(device::Device; dirichlet_boundary = [], periodic_coupling = [])

    function bilinear_kernel!(result, input, qpinfo)

        # the current cell region
        cell_region = qpinfo.region

        ε = tensor_view(input, 1, TDVector(6))
        σ = tensor_view(result, 1, TDVector(6))

        # simple linear elasticity problem
        mul!(σ, device.material_tensors[cell_region], ε)

        return nothing
    end

    # rhs resembles the various input strains
    function linear_kernel!(result, qpinfo)

        # the current cell region
        cell_region = qpinfo.region

        σ = tensor_view(result, 1, TDVector(6))

        # we assemble pre_stress data with a positive sign to the stress
        # since this is data in the linear kernel, it has to go to the rhs with a negative sign
        σ .= -device.pre_stress[cell_region]

        return nothing
    end

    PD = ProblemDescription("linear elastic problem")
    u = Unknown("displacement")
    assign_unknown!(PD, u)

    assign_operator!(PD, BilinearOperator(bilinear_kernel!, [εV(u, 1.0)], parallel = true))
    assign_operator!(PD, LinearOperator(linear_kernel!, [εV(u, 1.0)], parallel = true))

    # add dirichlet boundary data
    for db in dirichlet_boundary
        # assign_restriction!(
        #     PD,
        #     BoundaryDataRestriction(u; regions = [db.first], db.second)
        # )
        assign_operator!(
            PD,
            HomogeneousBoundaryData(u; regions = [db.first], value = db.second)
        )
    end

    # add periodic coupling
    for p in periodic_coupling
        assign_restriction!(PD, CoupledDofsRestriction(u, p..., parallel = true))
    end

    return PD
end
