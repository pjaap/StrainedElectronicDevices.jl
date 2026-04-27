struct SiO₂ <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    ε_r::Float64
    CTE::Float64
end


"""
    SiO₂()

    return the material SiO₂
"""
function SiO₂()
    # WARNING SiO₂ is amorphous and the isotropic simplification used here may be unsuitable
    # values taken from DOI 10.1063/1.4928320

    E = 73 # GPa ± 19
    ν = 0.17 # ± 0.1

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

    CTE = 6.5e-7 # https://www.azom.com/properties.aspx?ArticleID=1114

    return SiO₂(
        "SiO₂",
        matrix,
        3.9,
        CTE
    )
end

"""
    SiO2()

    returns SiO₂()
    ... for those who are afraid of Unicode 😄
"""
SiO2() = SiO₂()
