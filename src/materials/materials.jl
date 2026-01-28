abstract type AbstractMaterial end

"""
    NoMaterial

    Dummy for a non existent material for unused cell regions.
"""
struct NoMaterial <: AbstractMaterial
    # we need a dummy material tensor
    C::SMatrix{6, 6}
end
NoMaterial() = NoMaterial(@SArray zeros(6, 6))
const no_material = NoMaterial()

"""
    material_vector(n::Int)

    create a vector containing 'n' materials, all initialized with
    NoMaterial() and ready for setting up certain material regions.
"""
function material_vector(n::Int)
    result = Vector{AbstractMaterial}(undef, n)
    fill!(result, no_material)
    return result
end

"""
    Eν2lamé(E, ν)

    Return Lamé parameters λ, μ converted from E and ν (Young's modulus & Poisson ratio)
"""
function Eν2lamé(E, ν)
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))
    return λ, μ
end
