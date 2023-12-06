module TestPermutationTest

using Test
using CSV
using DataFrames
using StableRNGs
using NegativeControl
using TMLE
using Arrow
using Serialization

TESTDIR = joinpath(pkgdir(NegativeControl), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test permuted_estimand!" begin
    estimates = make_estimates()
    Ψ = estimates[1].estimand
    @test Ψ isa TMLE.StatisticalIATE
    # Treatment and Outcome
    permutation_variables = Set([Ψ.outcome, :RSID_103])
    Ψpermuted = NegativeControl.permuted_estimand!(permutation_variables, Ψ)
    @test Ψpermuted.outcome == NegativeControl.permuted_name(Ψ.outcome)
    @test keys(Ψpermuted.treatment_values) == (:RSID_103_permuted, :rs10043934)
    @test values(Ψpermuted.treatment_values) == values(Ψ.treatment_values)
    @test keys(Ψpermuted.treatment_confounders) == (:RSID_103_permuted, :rs10043934)
    @test values(Ψpermuted.treatment_confounders) == values(Ψ.treatment_confounders)
    @test Ψpermuted.outcome_extra_covariates == Ψ.outcome_extra_covariates
    # Treatment only
    permutation_variables = Set([:rs10043934])
    Ψpermuted = NegativeControl.permuted_estimand!(permutation_variables, Ψ)
    @test Ψpermuted.outcome == Ψ.outcome
    @test keys(Ψpermuted.treatment_values) == (:RSID_103, :rs10043934_permuted)
    @test values(Ψpermuted.treatment_values) == values(Ψ.treatment_values)
    @test keys(Ψpermuted.treatment_confounders) == (:RSID_103, :rs10043934_permuted)
    @test values(Ψpermuted.treatment_confounders) == values(Ψ.treatment_confounders)
    @test Ψpermuted.outcome_extra_covariates == Ψ.outcome_extra_covariates
    # Outcome only
    permutation_variables = Set([Ψ.outcome])
    Ψpermuted = NegativeControl.permuted_estimand!(permutation_variables, Ψ)
    @test Ψpermuted.outcome == NegativeControl.permuted_name(Ψ.outcome)
    @test Ψpermuted.treatment_values == Ψ.treatment_values
    @test Ψpermuted.treatment_confounders == Ψ.treatment_confounders
    @test Ψpermuted.outcome_extra_covariates == Ψ.outcome_extra_covariates
    # Composed Estimand
    Ψ = estimates[3].estimand
    @test Ψ isa ComposedEstimand
    outcome = Symbol("High light scatter reticulocyte percentage")
    permutation_variables = Set([:RSID_103, outcome])
    Ψpermuted = NegativeControl.permuted_estimand!(permutation_variables, Ψ)
    arg₁ = Ψpermuted.args[1]
    arg₂ = Ψpermuted.args[2]
    @test Ψpermuted.f == Ψ.f
    @test arg₁.outcome == arg₂.outcome == NegativeControl.permuted_name(outcome)
    @test keys(arg₁.treatment_values) == keys(arg₂.treatment_values) == (:RSID_103_permuted, :rs10043934)
    @test values(arg₁.treatment_values) == values(Ψ.args[1].treatment_values)
    @test values(arg₂.treatment_values) == values(Ψ.args[2].treatment_values)
    @test arg₂.outcome_extra_covariates == arg₁.outcome_extra_covariates == Ψ.args[1].outcome_extra_covariates
end


@testset "Test make_permutation_parameters" begin
    estimands = [Ψ.estimand for Ψ ∈ make_estimates()]
    expected_permuted_variables = Set([
        :rs117913124,
        :rs10043934,
        :RSID_104,
        Symbol("L50-L54 Urticaria and erythema"),
        Symbol("High light scatter reticulocyte percentage"),
        :RSID_103
    ])
    # Each estimand has 2 treatment variables and 1 outcome
    # Each estimand will thus lead to 3 new permutation estimands
    # A total of 12
    permutation_estimands, all_permuted_variables = NegativeControl.make_permutation_parameters(estimands; optimize=false, orders=(1,))
    @test length(permutation_estimands) == 12
    @test all(values(Ψ.treatment_values) == values(estimands[1].treatment_values) for Ψ ∈ permutation_estimands[1:3])
    @test all(values(Ψ.treatment_values) == values(estimands[2].treatment_values) for Ψ ∈ permutation_estimands[4:6])
    @test all(Ψ isa ComposedEstimand for Ψ ∈ permutation_estimands[7:9])
    @test all(values(Ψ.treatment_values) == values(estimands[4].treatment_values) for Ψ ∈ permutation_estimands[10:12])

    @test all_permuted_variables == expected_permuted_variables
    # For each of the 4 parameters there are:
    # - 3 x order 1 permutations
    # - 3 x order 2 permutations
    # - 1 order 3 permutation
    # = 7 parameters
    # Total: 7*4 = 28 permutation estimands
    permutation_estimands, all_permuted_variables = NegativeControl.make_permutation_parameters(estimands; optimize=false, orders=(1, 2, 3))
    @test length(permutation_estimands) == 28
    @test all(values(Ψ.treatment_values) == values(estimands[1].treatment_values) for Ψ ∈ permutation_estimands[1:7])
    @test all(values(Ψ.treatment_values) == values(estimands[2].treatment_values) for Ψ ∈ permutation_estimands[8:14])
    @test all(Ψ isa ComposedEstimand for Ψ ∈ permutation_estimands[15:21])
    @test all(values(Ψ.treatment_values) == values(estimands[4].treatment_values) for Ψ ∈ permutation_estimands[22:28])

    @test all_permuted_variables == expected_permuted_variables
end

@testset "Test run_permutation_test" begin
    make_fake_outputs()
    parsed_args = Dict(
        "dataset" => joinpath(TESTDIR, "data", "final.data.csv"),
        "results" => "tmle_output.hdf5",
        "outdir" => ".",
        "pval-threshold" => 1e-10,
        "verbosity" => 0,
        "limit" => nothing,
        "rng" => 123,
        "orders" => "1,2",
        "chunksize" => 5
    )
    generate_permutation_parameters_and_dataset(parsed_args)

    # Check permutation dataset file
    data = DataFrame(Arrow.Table(joinpath(parsed_args["outdir"], "permutation_dataset.arrow")))
    permuted_cols = filter(x -> endswith(x, "permuted"), names(data))
    @test Set(permuted_cols) == Set([
        "rs10043934_permuted",
        "High light scatter reticulocyte percentage_permuted",
        "RSID_103_permuted"
    ])
    
    # Only two estimates pass the threshold
    # Each estimate produces 6 new estimands
    # Total: 12 split into 3 files of chunksize 5
    batch_1 = deserialize("permutation_estimands_1.yaml").estimands
    @test length(batch_1) == 5
    batch_2 = deserialize("permutation_estimands_2.yaml").estimands
    @test length(batch_2) == 5
    batch_3 = deserialize("permutation_estimands_3.yaml").estimands
    @test length(batch_3) == 2

    # Clean
    clean()
    for file in readdir()
        if startswith(file, "permutation_")
            rm(file)
        end
    end
end

end

true