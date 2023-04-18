module TestPermutationTest

using Test
using CSV
using DataFrames
using StableRNGs
using NegativeControl

results_file = joinpath("data", "summary.csv")

@testset "Test misc functions" begin
    # Test parse_case_control
    @test NegativeControl.parse_case_control("2_&_1") == [2, 1]
    @test NegativeControl.parse_case_control("AC_&_CC") == ["AC", "CC"]
    @test NegativeControl.parse_case_control("AC_&_2") == ["AC", 2]

    # Test permuted_name
    @test NegativeControl.permuted_name("toto") == "toto_permuted"

    # Test split_string
    @test NegativeControl.split_string("AC_&_CC") == ["AC", "CC"]

    # Test make_permuted_col!
    data = DataFrame(
        A = [1, 2, 3],
        B = [1, 2, 3]
    )
    NegativeControl.make_permuted_col!(data, :A; rng=StableRNG(1234))
    @test data.A_permuted == [2, 3, 1]
    @test !hasproperty(data, :B_permuted)

    # Test make_parameter
    param_row = (CASE = "2_&_AC", CONTROL = "1_&_CC", CONFOUNDERS="PC1_&_PC2", COVARIATES="Age_&_Sex")
    param = NegativeControl.make_parameter(param_row, "Disease", ["rs12345", "rs56789_permuted"])
    expected_param = NegativeControl.IATE(
        target      = :Disease,
        treatment   = (rs12345=(case=2, control=1), rs56789_permuted=(case="AC", control="CC")),
        confounders = [:PC1, :PC2],
        covariates  = [:Age, :Sex]
    )
    @test param.target == expected_param.target
    @test param.confounders == expected_param.confounders
    @test param.covariates == expected_param.covariates
    @test param.treatment == expected_param.treatment

    # Test permute_treatments
    treatments = ["rs1234", "rs456", "sex"]
    comb = ["disease", "rs456"]
    target = "disease"
    new_treatments, new_target = NegativeControl.permutation_setting(comb, treatments, target)
    @test new_treatments == ["rs1234", "rs456_permuted", "sex"]
    @test new_target == "disease_permuted"
    comb = ["rs1234"]
    new_treatments, new_target = NegativeControl.permutation_setting(comb, treatments, target)
    @test new_treatments == ["rs1234_permuted", "rs456", "sex"]
    @test new_target == target
end

@testset "Test retrieve_significant_results" begin
    # The last is an ATE and filtered
    results = NegativeControl.retrieve_significant_results(results_file; threshold=:PVALUE => 1)
    @test size(results) == (8, 15)
    # With another P-value column constraint
    results = NegativeControl.retrieve_significant_results(results_file; threshold=:ADJUSTED_PVALUE => 0.9)
    @test size(results) == (7, 15)
end

@testset "Test make_permutation_parameters" begin
    results = NegativeControl.retrieve_significant_results(results_file; threshold=:PVALUE => 1)
    # Only looking at first 4 rows corresponding to two different treatments 
    # and 3 targets
    # For each of the 4 parameters there are:
    # - 3 x order 1 permutations
    # - 3 x order 2 permutations
    # - 1 order 3 permutation
    # = 7 parameters
    # Total: 7*4 = 28 parameters
    # Those parameters should be generated in optimal estimation order 
    # to save computations
    tocheck = 4
    parameters = NegativeControl.make_permutation_parameters(results[1:tocheck, :])
    @test size(parameters, 1) == 7*4

    expected_order = [
    ((:rs10043934, :rs17216707), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934, :rs17216707), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934, :rs17216707), Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted")),
    ((:rs10043934, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage")),
    ((:rs10043934, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage")),
    ((:rs10043934, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934, :rs17216707_permuted), Symbol("O68 Labour and delivery complicated by foetal stress [distress]")),
    ((:rs10043934, :rs17216707_permuted), Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted")),
    ((:rs10043934_permuted, :rs17216707), Symbol("High light scatter reticulocyte percentage")),
    ((:rs10043934_permuted, :rs17216707), Symbol("High light scatter reticulocyte percentage")),
    ((:rs10043934_permuted, :rs17216707), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934_permuted, :rs17216707), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934_permuted, :rs17216707), Symbol("O68 Labour and delivery complicated by foetal stress [distress]")),
    ((:rs10043934_permuted, :rs17216707), Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted")),
    ((:rs10043934_permuted, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage")),
    ((:rs10043934_permuted, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage")),
    ((:rs10043934_permuted, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934_permuted, :rs17216707_permuted), Symbol("High light scatter reticulocyte percentage_permuted")),
    ((:rs10043934_permuted, :rs17216707_permuted), Symbol("O68 Labour and delivery complicated by foetal stress [distress]")),
    ((:rs10043934_permuted, :rs17216707_permuted), Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted")),
    ((:rs10420720, :rs10741657), Symbol("E20-E35 Disorders of other endocrine glands_permuted")),
    ((:rs10420720, :rs10741657_permuted), Symbol("E20-E35 Disorders of other endocrine glands")),
    ((:rs10420720, :rs10741657_permuted), Symbol("E20-E35 Disorders of other endocrine glands_permuted")),
    ((:rs10420720_permuted, :rs10741657), Symbol("E20-E35 Disorders of other endocrine glands")),
    ((:rs10420720_permuted, :rs10741657), Symbol("E20-E35 Disorders of other endocrine glands_permuted")),
    ((:rs10420720_permuted, :rs10741657_permuted), Symbol("E20-E35 Disorders of other endocrine glands")),
    ((:rs10420720_permuted, :rs10741657_permuted), Symbol("E20-E35 Disorders of other endocrine glands_permuted")),
    ]
    for index in eachindex(parameters)
        param = parameters[index]
        treatment_target = tuple(keys(param.treatment), param.target)
        @test treatment_target == expected_order[index]
        @test param.covariates == [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]
        @test param.confounders == [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
    end

end

@testset "Test run_permutation_test" begin
    parsed_args = Dict(
        "resultdir" => "data",
        "outdir" => joinpath("results", "negative_control"),
        "pval-col" => "PVALUE",
        "pval-threshold" => 0.05,
        "estimator-file" => joinpath("data", "estimator.yaml"),
        "verbosity" => 0,
        "limit" => 10,
        "rng" => 123
    )    
    NegativeControl.run_permutation_test(parsed_args)

    summary = CSV.read(joinpath(parsed_args["outdir"], "summary.csv"), DataFrame)
    @test size(summary) == (10, 19)

    @test summary.TREATMENTS == [
        "rs10043934_&_rs17216707",
        "rs10043934_&_rs17216707",
        "rs10043934_&_rs17216707",
        "rs10043934_&_rs17216707_permuted",
        "rs10043934_&_rs17216707_permuted",
        "rs10043934_&_rs17216707_permuted",
        "rs10043934_&_rs17216707_permuted",
        "rs10043934_&_rs17216707_permuted",
        "rs10043934_&_rs17216707_permuted",
        "rs10043934_permuted_&_rs17216707"]
    @test summary.TARGET == [
        "High light scatter reticulocyte percentage_permuted",
        "High light scatter reticulocyte percentage_permuted",
        "O68 Labour and delivery complicated by foetal stress [distress]_permuted",
        "High light scatter reticulocyte percentage",
        "High light scatter reticulocyte percentage",
        "High light scatter reticulocyte percentage_permuted",
        "High light scatter reticulocyte percentage_permuted",
        "O68 Labour and delivery complicated by foetal stress [distress]",
        "O68 Labour and delivery complicated by foetal stress [distress]_permuted",
        "High light scatter reticulocyte percentage"
    ] 
    # All TMLEs have succeeded
    @test summary.TMLE_ESTIMATE isa Vector{Float64}
    @test all(x === missing for x in summary.LOG)
    # Clean
    rm(parsed_args["outdir"], force=true, recursive=true)
end

end

true