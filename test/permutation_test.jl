module TestPermutationTest

using Test
using CSV
using DataFrames
using StableRNGs
using NegativeControl
using TMLE
using Arrow
using Serialization

results_file = joinpath("data", "summary.csv")

@testset "Test misc functions" begin
    # Test parse_case_control
    @test NegativeControl.parse_case_control("2_&_1") == [2, 1]
    @test NegativeControl.parse_case_control("AC_&_CC") == ["AC", "CC"]
    @test NegativeControl.parse_case_control("AC_&_2") == ["AC", 2]

    # Test permuted_name
    @test NegativeControl.permuted_name("toto") == "toto_permuted"

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


@testset "Test make_permutation_parameters" begin
    results = NegativeControl.retrieve_significant_results(results_file; threshold=:PVALUE => 1)
    # Only looking at first 4 rows corresponding to two different treatments 
    # and 3 targets
    tocheck = 4
    results = results[1:tocheck, :]
    # First check for order 1 combinations: 
    # For each parameter there are 3 possible permutations
    # = 12 paramters
    parameters = NegativeControl.make_permutation_parameters(results, [1])
    expected_parameters = [
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934_permuted = (case = 2, control = 0), rs17216707 = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934_permuted = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]"), (rs10043934_permuted = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934 = (case = 2, control = 0), rs17216707 = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934 = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted"), (rs10043934 = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934 = (case = 2, control = 0), rs17216707_permuted = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934 = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]"), (rs10043934 = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands"), (rs10420720_permuted = (case = 1, control = 0), rs10741657 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands_permuted"), (rs10420720 = (case = 1, control = 0), rs10741657 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands"), (rs10420720 = (case = 1, control = 0), rs10741657_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
    ]
    @test parameters == expected_parameters

    # For each of the 4 parameters there are:
    # - 3 x order 1 permutations
    # - 3 x order 2 permutations
    # - 1 order 3 permutation
    # = 7 parameters
    # Total: 7*4 = 28 parameters
    # Those parameters should be generated in optimal estimation order 
    # to save computations
    parameters = NegativeControl.make_permutation_parameters(results, [1,2,3])
    @test size(parameters, 1) == 7*4
    expected_parameters = [
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934_permuted = (case = 2, control = 0), rs17216707 = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934_permuted = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934_permuted = (case = 2, control = 0), rs17216707 = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934_permuted = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]"), (rs10043934_permuted = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted"), (rs10043934_permuted = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934_permuted = (case = 2, control = 0), rs17216707_permuted = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934_permuted = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934_permuted = (case = 2, control = 0), rs17216707_permuted = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934_permuted = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]"), (rs10043934_permuted = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted"), (rs10043934_permuted = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934 = (case = 2, control = 0), rs17216707 = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934 = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted"), (rs10043934 = (case = 1, control = 0), rs17216707 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934 = (case = 2, control = 0), rs17216707_permuted = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage"), (rs10043934 = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934 = (case = 2, control = 0), rs17216707_permuted = (case = 1, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("High light scatter reticulocyte percentage_permuted"), (rs10043934 = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]"), (rs10043934 = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("O68 Labour and delivery complicated by foetal stress [distress]_permuted"), (rs10043934 = (case = 1, control = 0), rs17216707_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands"), (rs10420720_permuted = (case = 1, control = 0), rs10741657 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands_permuted"), (rs10420720_permuted = (case = 1, control = 0), rs10741657 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands"), (rs10420720_permuted = (case = 1, control = 0), rs10741657_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands_permuted"), (rs10420720_permuted = (case = 1, control = 0), rs10741657_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands_permuted"), (rs10420720 = (case = 1, control = 0), rs10741657 = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands"), (rs10420720 = (case = 1, control = 0), rs10741657_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
        IATE(Symbol("E20-E35 Disorders of other endocrine glands_permuted"), (rs10420720 = (case = 1, control = 0), rs10741657_permuted = (case = 2, control = 0)), [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6], [Symbol("Age-Assessment"), Symbol("Genetic-Sex")]),
    ]
    @test parameters == expected_parameters
end

@testset "Test run_permutation_test" begin
    outdir = "resultdir"
    ispath(outdir) || mkdir(outdir)
    parsed_args = Dict(
        "dataset" => joinpath("data", "final.data.csv"),
        "results" => joinpath("data", "summary.csv"),
        "outdir" => outdir,
        "pval-col" => "PVALUE",
        "pval-threshold" => 0.05,
        "verbosity" => 0,
        "limit" => 12,
        "rng" => 123,
        "orders" => "1",
        "chunksize" => 5
    )    
    generate_permutation_parameters_and_dataset(parsed_args)

    # Check permutation dataset file
    data = DataFrame(Arrow.Table(joinpath(parsed_args["outdir"], "permutation_dataset.arrow")))
    source_data = CSV.read(parsed_args["dataset"], DataFrame)
    results = NegativeControl.retrieve_significant_results(
        parsed_args["results"]; 
        threshold=Symbol(parsed_args["pval-col"]) => parsed_args["pval-threshold"]
    )
    treatment_cols = unique(reduce(vcat, [split(x, "_&_") for x in results[!, :TREATMENTS]]))
    target_cols = unique(results.TARGET)

    for colname in vcat(treatment_cols, target_cols)
        @test data[!, colname] !== data[!, NegativeControl.permuted_name(colname)]
    end
    
    # Check parameter files
    all_parameters = TMLE.Parameter[]
    for index in 1:3
        params_group = deserialize(joinpath(outdir, string("permutation_param_", index, ".bin")))
        append!(all_parameters, params_group)
    end
    @test length(all_parameters) == 12
    # Clean
    rm(outdir, force=true, recursive=true)
end

end

true