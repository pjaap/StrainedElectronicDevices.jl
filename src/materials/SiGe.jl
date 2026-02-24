# this file is based on pdelib2 (https://www.wias-berlin.de/software/index.jsp?id=pdelib)
# project "qubus"

struct SiGe <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    a::Float64
    c::Float64
    CTE::Float64
    ε_r::Float64
end


"""
    SiGe( x )

    return the material Si₁₋ₓGeₓ
"""
function SiGe(x)

    C11 = (165.8 - 37.3 * x)
    C12 = (63.9 - 15.6 * x)
    C44 = (79.6 - 12.8 * x)

    matrix = @SArray [
        C11 C12 C12 0   0   0
        C12 C11 C12 0   0   0
        C12 C12 C11 0   0   0
        0   0   0   C44 0   0
        0   0   0   0   C44 0
        0   0   0   0   0   C44
    ]

    CTE_Si = 0.76e-6
    CTE_Ge = 3.293e-6

    return SiGe(
        "SiGe",
        matrix,
        (5.431 + 0.2x + 0.027x^2),
        (5.431 + 0.2x + 0.027x^2),
        (1 - x) * CTE_Si + x * CTE_Ge,
        11.7 + 4.5x
    )
end

"""
    Si()

    Pure Silicon, alias for SiGe(0.0)
"""
Si() = SiGe(0.0)


"""
    Ge()

    Pure Germanium, alias for SiGe(1.0)
"""
Ge() = SiGe(1.0)
