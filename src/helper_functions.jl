"""
compute_lattice_mismatch(ref_material, material)

Compute the pre-strain induced by lattice mismatch between reference material and a given material.
The result is given in Voigt notation
"""
function compute_lattice_mismatch(ref_material, material)
    return @SArray [
        (ref_material.a - material.a) / material.a, # ε_xx
        (ref_material.a - material.a) / material.a, # ε_yy
        (ref_material.c - material.c) / material.c, # ε_zz
        0, 0, 0, # Voigt off diagonal
    ]
end


"""
    z_rotation_matrix(α)

    Rotation matrix (ℝ³ˣ³) for rotation around the z axis by angle α
"""
function z_rotation_matrix(α)
    return @SArray [
        cos(α) -sin(α) 0.0
        sin(α)  cos(α) 0.0
        0.0     0.0    1.0
    ]
end


"""
    y_rotation_matrix(α)

    Rotation matrix (ℝ³ˣ³) for rotation around the y axis by angle α
"""
function y_rotation_matrix(α)
    return @SArray [
        cos(α) 0.0 -sin(α)
        0.0    1.0  0.0
        sin(α) 0.0  cos(α)
    ]
end


## transform Voigt vector to matrix
function voigt2strain(εv)
    return @SArray [
        εv[1]     εv[6] / 2 εv[5] / 2
        εv[6] / 2 εv[2]     εv[4] / 2
        εv[5] / 2 εv[4] / 2 εv[3]
    ]
end
function voigt2stress(σv)
    return @SArray [
        σv[1] σv[6] σv[5]
        σv[6] σv[2] σv[4]
        σv[5] σv[4] σv[3]
    ]
end


## transform matrix to Voigt vector
strain2voigt(ε) = @SArray [ε[1, 1], ε[2, 2], ε[3, 3], 2ε[2, 3], 2ε[1, 3], 2ε[1, 2]]
stress2voigt(σ) = @SArray [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]

"""
    get_rotated_tensor(𝐂, R)

    Compute a material tensor representing a rotated crystal structure, as given by the rotation matrix R.
    This linear mapping is computed by AD.
"""
function get_rotated_tensor(𝐂ⱽ, R)

    # define the linear mapping in matrix representation in order to apply the rotation
    function mapping_eval(εV) # Voigt Strain
        return stress2voigt(R' * voigt2stress(𝐂ⱽ * strain2voigt(R * voigt2strain(εV) * R')) * R)
    end

    # compute the rotated linear mapping
    return SMatrix{6, 6}(ForwardDiff.jacobian(mapping_eval, zeros(6)))
end


"""
    rotate_crystal_structure!(material::AbstractMaterial, R::AbstractMatrix)

    Rotate the intrinsic crystal structure of the given material by the rotation matrix R.
    The rotated mapping of the material tensor 𝐂 is then
        σ = Rᵀ 𝐂[ R ε Rᵀ ] R
"""
function rotate_crystal_structure(material::T, R::AbstractMatrix) where {T <: AbstractMaterial}
    𝐂 = get_rotated_tensor(material.C, R)

    # What do the following lines mean?
    # We want to have immutable structs for the materials. period. So construct a new material "T"
    # However, each material has different field names.
    # Hence, we have to filter out the :C symbol and replace it with the rotated tensor.
    data = [ [field, getproperty(material, field)] for field in propertynames(material) ]

    # manipulate :C
    data[findfirst(d -> d[1] == :C, data)][2] = 𝐂

    # also modify the name, to make sure that the user knows that this is not the original material...
    data[findfirst(d -> d[1] == :name, data)][2] = "rotated_" * material.name

    # construct the type with manipulated data
    return T(last.(data)...)
end
