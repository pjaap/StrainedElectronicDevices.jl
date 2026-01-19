struct SiO₂ <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
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

    return SiO₂("SiO₂", matrix)
end

"""
    SiO2()

    returns SiO₂()
    ... for those who are afraid of Unicode 😄
"""
SiO2() = SiO₂()
