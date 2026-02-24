module StrainedElectronicDevices

using ExtendableFEM: ExtendableFEM, BilinearOperator, LinearOperator,
    ProblemDescription, TDVector, Unknown, assign_operator!,
    assign_unknown!, tensor_view, εV, assign_restriction!, BoundaryDataRestriction, HomogeneousBoundaryData,
    CoupledDofsRestriction, grad, HomogeneousData
using ExtendableGrids: ExtendableGrid, num_cellregions, dim_space
import ForwardDiff
using LinearAlgebra: mul!
using StaticArrays: @SArray, SMatrix, SVector
using SimplexGridFactory: SimplexGridFactory, SimplexGridBuilder, cellregion!,
    facet!, facetregion!, maxvolume!, point!, regionpoint!, simplexgrid
using TetGen: TetGen


gridsdir(args...) = joinpath(pkgdir(StrainedElectronicDevices), "src", "grids", args...)
scriptsdir(args...) = joinpath(pkgdir(StrainedElectronicDevices), "scripts", args...)


include("materials/materials.jl")
export AbstractMaterial, NoMaterial, material_vector, no_material

# materials
include("materials/Al.jl")
include("materials/Al2O3.jl")
include("materials/AlGaN.jl")
include("materials/Pd.jl")
include("materials/SiGe.jl")
include("materials/SiO2.jl")
include("materials/TiN.jl")
export AlGaN, GaN, AlN, Si, Ge, SiGe, SiO2, SiO₂, TiN, Al₂O₃, Al2O3, Al, Pd


include("grids/cuboid_with_QW.jl")
include("grids/cuboid_with_wire.jl")
include("grids/grid_Corley_Wiciak.jl")
export create_cuboid_grid_with_quantum_well
export create_cuboid_grid_with_wire
export grid_Corley_Wiciak_3D

include("device.jl")
export Device

include("linear_elasticity_problem.jl")
export create_linear_elasticity_problem

include("electrostatic_problem.jl")
export create_electrostatic_problem

include("helper_functions.jl")
export compute_lattice_mismatch
export z_rotation_matrix
export y_rotation_matrix
export rotate_crystal_structure


end # module
