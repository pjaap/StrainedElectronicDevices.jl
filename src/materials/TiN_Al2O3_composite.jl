# TiN / Al₂O₃ composite (50/50 vol. fraction) – transversally isotropic

struct TiN_Al2O3_composite <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    CTE::Float64
end


"""
    TiN_Al2O3_composite()

    Mixed layer of TiN (metal) and Al₂O₃ (insulator) with ~50/50 volume fraction.
    Transversally isotropic homogenized elasticity tensor.

    C (GPa):
        C₁₁ = 216.8, C₁₂ = 68.4, C₁₃ = 61.0, C₃₃ = 187.9,
        C₄₄ = 63.8, C₆₆ = 74.2
    α_CTE = 7.7 × 10⁻⁶ 1/K (average between room temp and cryogenic)
"""
function TiN_Al2O3_composite()

    C = @SArray [
        216.8  68.4  61.0   0     0    0
        68.4  216.8 61.0   0     0    0
        61.0  61.0  187.9   0     0    0
        0      0    0     63.8  0    0
        0      0    0      0    63.8 0
        0      0    0      0     0   74.2
    ]

    CTE = 7.7e-6  # 1/K

    return TiN_Al2O3_composite(
        "TiN/Al₂O₃ composite (50/50)",
        C,
        CTE
    )
end
