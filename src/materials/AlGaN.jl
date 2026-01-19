# this file is based on pdelib2 (https://www.wias-berlin.de/software/index.jsp?id=pdelib)
# project "qubus"

struct AlGaN <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    a::Float64
    c::Float64
    E::SMatrix{3, 6}
    ε_r::Float64
end

"""
    AlGaN(x)

    return the material Al₁₋ₓGaₓN with x ∈ [0,1]
"""
function AlGaN(x)

    aGaN = 3.189
    cGaN = 5.185
    c11GaN = 390
    c12GaN = 145
    c13GaN = 106
    c33GaN = 398
    c44GaN = 105
    e31GaN = -0.49
    e33GaN = 0.73
    e15GaN = 0
    pSpGaN = -0.034
    epsGaN = 9.8

    aAlN = 3.112
    cAlN = 4.982
    c11AlN = 396
    c12AlN = 137
    c13AlN = 108
    c33AlN = 373
    c44AlN = 116
    e31AlN = -0.44
    e33AlN = 0.744
    e15AlN = 0
    pSpAlN = -0.09
    epsAlN = 10.1

    aAlGaN = (1 - x) * aAlN + x * aGaN
    cAlGaN = (1 - x) * cAlN + x * cGaN
    c11 = (1 - x) * c11AlN + x * c11GaN
    c12 = (1 - x) * c12AlN + x * c12GaN
    c13 = (1 - x) * c13AlN + x * c13GaN
    c33 = (1 - x) * c33AlN + x * c33GaN
    c44 = (1 - x) * c44AlN + x * c44GaN
    e31 = (1 - x) * e31AlN + x * e31GaN
    e33 = (1 - x) * e33AlN + x * e33GaN
    e15 = (1 - x) * e15AlN + x * e15GaN

    ε_r = (1 - x) * epsAlN + x * epsGaN

    # stiffness tensor
    matrixC = @SArray [
        c11 c12 c13  0   0   0
        c12 c11 c13  0   0   0
        c13 c13 c33  0   0   0
        0   0   0    c44 0   0
        0   0   0    0   c44 0
        0   0   0    0   0   c44
    ]

    # piezoelectric tensor
    matrixE = @SArray [
        0   0   0   0   e15 0
        0   0   0   e15 0   0
        e31 e31 e33 0   0   0
    ]

    return AlGaN(
        "wz_Al₁₋ₓGaₓN",
        matrixC,
        aAlGaN,
        cAlGaN,
        matrixE,
        ε_r,
    )

    return material
end

"""
    GaN()

    convenient alias for AlGaN(1.0)
"""
GaN() = AlGaN(1.0)

"""
    AlN()

    convenient alias for AlGaN(0.0)
"""
AlN() = AlGaN(1.0)
