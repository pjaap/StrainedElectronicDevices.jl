# β-Si₃N₄ (beta silicon nitride) – hexagonal (transversally isotropic) elasticity

struct Si₃N₄ <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    CTE::Float64
end


"""
    Si₃N₄()

    β-Si₃N₄ (beta silicon nitride, hexagonal, transversally isotropic).

    Elasticity tensor (GPa):
        C11 = 407, C12 = 178, C13 = 104, C33 = 407,
        C44 = 100, C66 = 115
    α  = 0.5 × 10⁻⁶ 1/K
    Biaxial pre-stress σ* = diag(1, 1, 0) GPa (tensile)
    Source: Woods
"""
function Si₃N₄()

    C11 = 407.0  # GPa
    C12 = 178.0  # GPa
    C13 = 104.0  # GPa
    C33 = 407.0  # GPa
    C44 = 100.0  # GPa
    C66 = 115.0  # GPa

    matrix = @SArray [
        C11 C12 C13  0    0    0
        C12 C11 C13  0    0    0
        C13 C13 C33  0    0    0
        0   0   0   C44   0    0
        0   0   0    0   C44   0
        0   0   0    0    0   C66
    ]

    CTE = 0.5e-6  # 1/K (average thermal expansion coefficient)

    return Si₃N₄(
        "β-Si₃N₄ (hexagonal)",
        matrix,
        CTE
    )
end

"""
    Si3N4()

    returns Si₃N₄() — for those who are afraid of Unicode 😄
"""
Si3N4() = Si₃N₄()
