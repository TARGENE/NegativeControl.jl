module TestPermutationTest

using Test
using CSV
using DataFrames
using StableRNGs

@testset "Test misc functions" begin
    # Test parse_case_control
    @test parse_case_control("2_&_1") == [2, 1]
    @test parse_case_control("AC_&_CC") == ["AC", "CC"]
    @test parse_case_control("AC_&_2") == ["AC", 2]

    # Test permuted_name
    @test permuted_name("toto") == "toto_permuted"

    # Test split_string
    @test split_string("AC_&_CC") == ["AC", "CC"]

    # Test make_permuted_col!
    data = DataFrame(
        A = [1, 2, 3],
        B = [1, 2, 3]
    )
    make_permuted_col!(data, :A; rng=StableRNG(1234))
    @test data.A_permuted == [2, 3, 1]
    @test !hasproperty(data, :B_permuted)

    # Test make_parameter
    param_row = (CASE = "2_&_AC", CONTROL = "1_&_CC", CONFOUNDERS="PC1_&_PC2", COVARIATES="Age_&_Sex")
    param = make_parameter(param_row, "Disease", ["rs12345", "rs56789_permuted"])
    expected_param = IATE(
        target      = :Disease,
        treatment   = (rs12345=(case=2, control=1), rs56789_permuted=(case="AC", control="CC")),
        confounders = [:PC1, :PC2],
        covariates  = [:Age, :Sex]
    )
    @test param.target == expected_param.target
    @test param.confounders == expected_param.confounders
    @test param.covariates == expected_param.covariates
    @test param.treatment == expected_param.treatment
    @test param isa IATE

    # Test permute_treatments
    treatments = ["rs1234", "rs456", "sex"]
    treatment_comb = ["sex", "rs456"]
    permuted_treatments = permute_treatments(treatments, treatment_comb)
    @test permuted_treatments == ["rs1234", "rs456_permuted", "sex_permuted"]
    treatment_comb = ["rs1234"]
    permuted_treatments = permute_treatments(treatments, treatment_comb)
    @test permuted_treatments == ["rs1234_permuted", "rs456", "sex"]
end

@testset "Test retrieve_significant_results" begin
    results_file = joinpath("test", "data", "summary.csv")
    # The last is an ATE and filtered
    results = retrieve_significant_results(joinpath("test", "data", "summary.csv"); threshold=:PVALUE => 1)
    @test size(results) == (8, 15)
    # With another P-value column constraint
    results = retrieve_significant_results(joinpath("test", "data", "summary.csv"); threshold=:ADJUSTED_PVALUE => 0.9)
    @test size(results) == (7, 15)
end

@testset "Test make_permutation_parameters" begin
    results = retrieve_significant_results(joinpath("test", "data", "summary.csv"); threshold=:PVALUE => 1)
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
    parameters = make_permutation_parameters(results[1:tocheck, :])
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

@testset "Test run" begin
    parsed_args = Dict(
        "resultdir" => "data",
        "outdir" => pwd(),
        "pval-col" => "PVALUE",
        "pval-threshold" => 0.05,
        "estimator-file" => joinpath("data", "estimator.yaml"),
        "verbosity" => 1,
        "limit" => 10
    )    


    parsed_args = Dict(
        "resultdir" => "results",
        "outdir" => "results/negative_control",
        "pval-col" => "PVALUE",
        "pval-threshold" => 0.05,
        "estimator-file" => joinpath("data", "estimator.yaml"),
        "verbosity" => 1,
        "limit" => 10
    )    
    # run(parsed_args)
end

end

true