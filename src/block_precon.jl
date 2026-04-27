struct RestrictedBlockPreconditioner{FT}
    restricted_indices::AbstractVector
    partitions_free_indices::Vector # length N
    factorization_restricted::FT
    factorizations_free::Vector{FT} # length N
end


function RestrictedBlockPreconditioner(M::AbstractMatrix; blocks = 1, factorization = LinearAlgebra.lu, stabilizer = 0.0, verbosity = 0)

    # we have M = [ A  B ]
    #             [ Bᵀ 0 ]

    # we assume that A has nonzero diagonal
    dofs_first_block = findlast(≠(0), diag(M))

    # I see no other way than creating this expensive copy
    B = M[1:dofs_first_block, (dofs_first_block + 1):end]

    # find all rows with entries in B ⇒ these are the restricted indices
    restricted_indices, _, _ = findnz(B)

    # ... plus the col indices of B!
    restricted_indices = vcat(unique(restricted_indices), (dofs_first_block + 1):size(M, 1))

    free_indices = setdiff(1:size(M, 2), restricted_indices)

    partitions_free_indices = collect(chunks(free_indices, n = blocks))

    M_restricted = M[restricted_indices, restricted_indices]
    M_restricted_fac = factorization(M_restricted)

    factorizations_free = typeof(M_restricted_fac)[]
    for part in partitions_free_indices
        M_part_fac = factorization(M[part, part] + stabilizer * LinearAlgebra.I)
        push!(factorizations_free, M_part_fac)
    end

    return RestrictedBlockPreconditioner(restricted_indices, partitions_free_indices, M_restricted_fac, factorizations_free)
end


function RestrictedBlockPreconBuilder(factorization = LinearAlgebra.lu; blocks = 1, verbosity = 0, stabilizer = 0.0)

    # this is the resulting LinearSolve.jl compatible preconditioner
    function prec(M, p)
        RBPC = RestrictedBlockPreconditioner(M; blocks, factorization, stabilizer, verbosity)
        return (RBPC, LinearAlgebra.I)
    end

    return prec
end


function LinearAlgebra.ldiv!(u, p::RestrictedBlockPreconditioner, v)

    t1 = Threads.@spawn  begin
        # sigh, lu cannot work on views. hard copy the data
        # TODO resolve this problem!
        u_part = u[p.restricted_indices]
        v_part = v[p.restricted_indices]
        @views LinearAlgebra.ldiv!(u_part, p.factorization_restricted, v_part)
        u[p.restricted_indices] = u_part
        v[p.restricted_indices] = v_part
    end

    Threads.@threads for i in eachindex(p.partitions_free_indices)
        part, fac = p.partitions_free_indices[i], p.factorizations_free[i]
        u_part = u[part]
        v_part = v[part]
        @views LinearAlgebra.ldiv!(u_part, fac, v_part)
        u[part] = u_part
        v[part] = v_part
    end

    fetch(t1)

    return u
end
