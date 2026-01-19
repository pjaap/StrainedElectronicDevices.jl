# this file is based on pdelib2 (https://www.wias-berlin.de/software/index.jsp?id=pdelib)
# project "qubus"

struct Al₂0₃ <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    CTE::Float64
end

"""
    Al₂0₃(x)

    return the material Al₂0₃
"""
function Al₂O₃()

    # cf. https://web.archive.org/web/20131029202129/http://www.ceramtec.de/files/ma_materials_data_de_en.pdf
    E = 390 # GPa
    ν = 0.23

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

    return Al₂0₃(
        "Al₂0₃",
        matrix,
        3.3e-6 # CTE
    )
end

"""
    Al2O3()

    returns Al₂0₃()
    ... for those who are afraid of Unicode 😄
"""
Al2O3() = Al₂0₃()
