using Aqua
using ExplicitImports
using StrainedElectronicDevices
using Test

@testset "Quality Tests" begin
    @testset "Aqua.jl" begin
        Aqua.test_all(StrainedElectronicDevices)
    end

    @testset "ExplicitImports.jl" begin
        @test ExplicitImports.check_all_explicit_imports_are_public(StrainedElectronicDevices) === nothing
        @test ExplicitImports.check_all_explicit_imports_via_owners(StrainedElectronicDevices) === nothing
        @test ExplicitImports.check_all_qualified_accesses_are_public(
            StrainedElectronicDevices,
            ignore = (:jacobian,)
        ) === nothing
        @test ExplicitImports.check_all_qualified_accesses_via_owners(StrainedElectronicDevices) === nothing
        @test ExplicitImports.check_no_implicit_imports(StrainedElectronicDevices) === nothing
        @test ExplicitImports.check_no_self_qualified_accesses(StrainedElectronicDevices) === nothing
        @test ExplicitImports.check_no_stale_explicit_imports(StrainedElectronicDevices) === nothing
    end
end
