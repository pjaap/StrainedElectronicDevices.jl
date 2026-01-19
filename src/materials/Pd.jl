# this file is based on pdelib2 (https://www.wias-berlin.de/software/index.jsp?id=pdelib)
# project "qubus"

struct Pd <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    CTE::Float64
end

"""
    Pd(x)

    return the material Pd
"""
function Pd()

    # from https://en.wikipedia.org/wiki/Palladium
    E = 121 # GPa
    ν = 0.39

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
        "Pd",
        matrix,
        8.86e-6 # CTE
    )
end
