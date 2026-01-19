# this file is based on pdelib2 (https://www.wias-berlin.de/software/index.jsp?id=pdelib)
# project "qubus"

struct Al <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    CTE::Float64
end

"""
    Al(x)

    return the material Al
"""
function Al()

    # from https://en.wikipedia.org/wiki/Aluminium
    E = 70.0 # GPa
    ν = 0.35

    # Lamé parameters
    λ, μ = Eν2lamé(E, ν)

    C11 = λ + 2μ
    C12 = λ
    C44 = μ

    matrix = @SArray [
        C11 C12 C12 0   0   0
        C12 C11 C12 0   0   0
        C12 C12 C11 0   0   0
        0   0   0   C44 0   0
        0   0   0   0   C44 0
        0   0   0   0   0   C44
    ]

    return Al(
        "Al",
        matrix,
        14.16e-6 # CTE
    )
end
