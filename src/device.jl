struct Device{Tc, Ti}

    grid::ExtendableGrid{Tc, Ti}

    pre_stress::Vector{SVector{6}}

    material_tensors::Vector{SMatrix{6, 6}}

    dielectric_permittivity::Vector{Float64}
end


"""
    Device(
        grid::ExtendableGrid,
        materials::Vector{AbstractMaterial};
        thermal_strain = [],
        pre_stress = [],
        pre_strain = [],
        lattice_mismatch = []
    )

    Create a Device with specified geometry, material layers and physical properties.

    grid: ExtendableGrids compatible grid. Materials are given by cell regions
    materials: Vector of AbstractMaterial, according to the grid cell regions

    optional parameters:
      thermal_strain: Vector{Pair} defining thermal strains for grid cell regions in Voigt notation
      pre_stress: Vector{Pair} defining given pre stresses for grid cell regions in Voigt notation
      pre_strain: Vector{Pair} defining thermal pre strains for grid cell regions in Voigt notation
      lattice_mismatch: Vector{Pair} defining lattice mismatches for grid cell regions in Voigt notation

    The optional parameters are given in the form , e.g, [ 2 => [ 1, 1, 1, 0, 0, 0], 4 => [ 1, 2, 3, 4, 5, 6] ]

    Note: if both pre_strain and pre_stress are given on the same cell region, only the pre_stress is considered.
"""
function Device(
        grid::ExtendableGrid{Tc, Ti},
        materials::Vector{<:AbstractMaterial};
        thermal_strain = [],
        pre_stress = [],
        pre_strain = [],
        lattice_mismatch = [],
        grid_scaling = 1.0
    ) where {Tc, Ti}

    cell_regions = num_cellregions(grid)
    dim = dim_space(grid)

    @assert dim == 3
    dim_Voigt = 6

    # prepare physical properties for each cell region
    ps = [ @SArray zeros(dim_Voigt) for i in 1:cell_regions ]
    mt = [ material.C for material in materials ]

    for (k, v) in pre_strain
        # convert strain to stress via material tensor as we are only interested in stress values
        ps[k] += mt[k] * v
    end

    for (k, v) in pre_stress
        ps[k] += v
    end

    for (k, v) in thermal_strain
        # convert strain to stress via material tensor as we are only interested in stress values
        ps[k] += mt[k] * v
    end

    for (k, v) in lattice_mismatch
        # convert strain to stress via material tensor as we are only interested in stress values
        ps[k] += mt[k] * v
    end

    ε_0	= 8.854187817620389e-12
    dielectric_permittivity = [ material.ε_r * ε_0 * grid_scaling for material in materials ]

    return Device{Tc, Ti}(grid, ps, mt, dielectric_permittivity)
end
