function create_electrostatic_problem(device::Device; dirichlet_regions = [], periodic_coupling = [], charge_density = x -> 0.0)

    function bilinear_kernel!(result, input, qpinfo)

        # the current cell region
        cell_region = qpinfo.region

        @. result = device.dielectric_permittivity[cell_region] * input

        return nothing
    end

    function linear_kernel!(result, qpinfo)

        result[1] = charge_density(qpinfo.x)

        return nothing
    end

    PD = ProblemDescription("electrostatic problem")
    ϕ = Unknown("electrostatic potential")
    assign_unknown!(PD, ϕ)

    assign_operator!(PD, BilinearOperator(bilinear_kernel!, [grad(ϕ)], parallel = true))
    assign_operator!(PD, LinearOperator(linear_kernel!, [grad(ϕ)], parallel = true))

    # add dirichlet boundary data
    for db in dirichlet_regions
        assign_operator!(
            PD,
            HomogeneousData(ϕ; regions = [db.first], value = db.second)
        )
    end

    # add periodic coupling
    for p in periodic_coupling
        assign_restriction!(PD, CoupledDofsRestriction(ϕ, p..., parallel = true))
    end

    return PD
end
