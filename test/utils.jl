module TestUtils

using Test
using NegativeControl

TESTDIR = joinpath(pkgdir(NegativeControl), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test treatment_variables and outcome_variables" begin
    estimates = make_estimates()
    # Classic estimand
    Ψ = estimates[1].estimand
    @test Ψ isa TMLE.StatisticalIATE
    @test NegativeControl.treatment_variables(Ψ) == [:RSID_103, :rs10043934]
    @test NegativeControl.outcome_variables(Ψ) == [Ψ.outcome]
    # Composed estimand
    Ψ = estimates[3].estimand
    @test Ψ isa ComposedEstimand
    @test NegativeControl.treatment_variables(Ψ) == [:RSID_103, :rs10043934]
    @test NegativeControl.outcome_variables(Ψ) == [Symbol("High light scatter reticulocyte percentage")]
end

@testset "Test read_significant_results" begin
    prefix = "tmle_output"
    estimates = make_estimates()
    save(estimates;prefix=prefix)

    threshold = 1.
    for ext in ("jls", "json", "hdf5")
        results = NegativeControl.read_significant_results(string(prefix, "." ,ext); threshold=threshold)
        @test length(results) == length(estimates)
        for (index, Ψ̋) ∈ enumerate(estimates)
            results[index] == Ψ̋.estimand
        end
    end

    threshold = 1e-9
    for ext in ("jls", "json", "hdf5")
        results = NegativeControl.read_significant_results(string(prefix, "." ,ext); threshold=threshold)
        @test length(results) == 2
    end

    threshold = 1e-100
    for ext in ("jls", "json", "hdf5")
        results = NegativeControl.read_significant_results(string(prefix, "." ,ext); threshold=threshold)
        @test length(results) == 0
    end

    clean()
end

end

true