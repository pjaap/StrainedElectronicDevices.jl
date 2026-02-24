struct TiN <: AbstractMaterial
    name::String
    C::SMatrix{6, 6}
    CTE::Float64
    ε_r::Float64 # this is only a dummy
end


"""
    TiN( model )

    model ∈ { :A, :B }
    The different models are explained in DOI 10.1103/PhysRevApplied.20.024056

    return the material TiN
"""
function TiN(model)

    # values from 10.1103/PhysRevApplied.20.024056
    if model == :A
        E = 106 # GPa ± 25
    elseif model == :B
        E = 91 # GPa ± 19
    end

    ν = 0.3 # ± 0.1

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

    CTE = 9.4 * 1.0e-6 # https://www.msesupplies.com/en-de/pages/list-of-thermal-expansion-coefficients-cte-for-natural-and-engineered-materials

    return TiN(
        "TiN Model $model",
        matrix,
        CTE,
        1.0
    )


end
