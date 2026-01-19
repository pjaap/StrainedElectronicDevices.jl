"""
    create_linear_elasticity_problem(device::Device, unknown::Unknown; kwargs...)

    This creates an ExtendableFEM compatible linear elasticity problem based on a given Device.

    TODO: specify the precise PDE
"""
function create_linear_elasticity_problem(device::Device, unknown::Unknown; kwargs...)

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

    PD = ProblemDescription()
    assign_unknown!(PD, unknown)

    assign_operator!(PD, BilinearOperator(bilinear_kernel!, [εV(unknown, 1.0)]))
    assign_operator!(PD, LinearOperator(linear_kernel!, [εV(unknown, 1.0)]))

    return PD
end
