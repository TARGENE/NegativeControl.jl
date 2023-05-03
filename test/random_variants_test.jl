module RandomVariantsTests

using Test
using NegativeControl
using BGEN
using Serialization
using StableRNGs
using HTTP
using Statistics
using TMLE

results_file = joinpath("data", "summary.csv")

@testset "Test same_maf" begin
    b = Bgen(
        BGEN.datadir("example.8bits.bgen"); 
        sample_path=BGEN.datadir("example.sample"), 
        idx_path=BGEN.datadir("example.8bits.bgen.bgi")
    )
    v = variant_by_rsid(b, "RSID_198")
    @test mean(minor_allele_dosage!(b, v)) ≈ 0.48411763 atol=1e-5
    @test NegativeControl.same_maf(b, v, 0.51; reltol=0.05) == false
    @test NegativeControl.same_maf(b, v, 0.50; reltol=0.05) == true
end

@testset "Test control_case" begin
    b = Bgen(
        BGEN.datadir("example.8bits.bgen"); 
        sample_path=BGEN.datadir("example.sample"), 
        idx_path=BGEN.datadir("example.8bits.bgen.bgi")
    )
    v₁ = variant_by_rsid(b, "RSID_198")
    minor_allele_dosage!(b, v₁)
    v₂ = variant_by_rsid(b, "RSID_199")
    minor_allele_dosage!(b, v₂)

    # Int representation
    @test NegativeControl.control_case(v₁, v₂, 0, 1) == (control=1, case=0)
    @test NegativeControl.control_case(v₁, v₂, 0, 2) == (control=2, case=0)
    @test NegativeControl.control_case(v₁, v₂, 2, 1) == (control=1, case=2)

    # String representation
    @test minor_allele(v₁) == "A"
    @test major_allele(v₁) == "G"
    @test minor_allele(v₂) == "G"
    @test major_allele(v₂) == "A"
    @test NegativeControl.control_case(v₁, v₂, "AG", "AA") == (control="GG", case="GA")
    @test NegativeControl.control_case(v₁, v₂, "GG", "AA") == (control="GG", case="AA")
    @test NegativeControl.control_case(v₁, v₂, "AA", "AG") == (control="GA", case="GG")
end

@testset "Test is_in_regulatory_region" begin
    v = deserialize(joinpath("data", "variants", "not_annotated_variant.bin"))
    @test NegativeControl.is_in_regulatory_region(v) == false

    v = deserialize(joinpath("data", "variants", "regulatory_variant.bin"))
    @test NegativeControl.is_in_regulatory_region(v) == true

    v = deserialize(joinpath("data", "variants", "failing_variant.bin"))
    @test_throws HTTP.Exceptions.StatusError NegativeControl.is_in_regulatory_region(v)
end

@testset "Test find_maf_matching_random_variants" begin
    # One transactor
    b = Bgen(
        BGEN.datadir("example.8bits.bgen"); 
        sample_path=BGEN.datadir("example.sample"), 
        idx_path=BGEN.datadir("example.8bits.bgen.bgi")
        )
    variant = variant_by_rsid(b, "RSID_103")
    transactors = Set(["RSID_103"])
    variants = NegativeControl.find_maf_matching_random_variants!(
        b, variant, transactors; 
        p=10, rng=StableRNG(123), reltol=0.05
    )
    @test size(variants, 1) == 10
    for v in variants
        @test rsid(v) != "RSID_103"
    end

    # All transactors
    p = 10
    bgen_prefix = joinpath("data", "bgen", "ukb")
    trans_actors = Set(["RSID_103", "RSID_198"])
    reltol = 0.1
    variant_map = NegativeControl.find_maf_matching_random_variants(
        trans_actors, 
        bgen_prefix,
        trans_actors; 
        p=p, rng=StableRNG(123), reltol=reltol
    )

    # Check criteria
    expected_mapped_rsids = Dict(
        "RSID_103" => [
            "RSID_110", "RSID_10", "RSID_79", "RSID_6", "RSID_179", 
            "RSID_142", "RSID_42", "RSID_106", "RSID_61", "RSID_57"
            ], 
        "RSID_198" => [
            "RSID_20", "RSID_58", "RSID_197", "RSID_120", "RSID_65", 
            "RSID_127", "RSID_98", "RSID_188", "RSID_165", "RSID_158",
            ]
    )
    for (key, mapped_rsids) in expected_mapped_rsids
        variant = variant_map[key][1]
        mapped_variants = variant_map[key][2]
        @test rsid(variant) == key
        # Ensures the process is deterministic
        @test [rsid(v) for v in mapped_variants] == mapped_rsids
        target_maf = mean(minor_allele_dosage!(b, variant))
        @test all(NegativeControl.same_maf(b, v, target_maf, reltol=reltol) for v in mapped_variants)
        @test all(!NegativeControl.is_in_regulatory_region(v) for v in mapped_variants)
    end

    # Failing 
    trans_actors = Set(["RSID_3"])
    reltol = 0.005
    @test_throws NegativeControl.NotEnoughMatchingVariantsError("RSID_3", p, reltol) NegativeControl.find_maf_matching_random_variants(
        trans_actors, 
        bgen_prefix,
        trans_actors; 
        p=p, rng=StableRNG(123), reltol=reltol
    )
end

@testset "Test make_random_variants_parameters" begin
    p = 5
    reltol = 0.1
    rng = StableRNG(123)
    bgen_prefix = joinpath("data", "bgen", "ukb")
    trans_actors = Set(["RSID_103", "RSID_198"])
    
    results = NegativeControl.retrieve_significant_results(results_file, threshold=:PVALUE => 0.05)
    variant_map = NegativeControl.find_maf_matching_random_variants(
        trans_actors, 
        bgen_prefix,
        trans_actors; 
        p=p, rng=rng, reltol=reltol
    )

    parameters = NegativeControl.make_random_variants_parameters(results, variant_map)
    # Only 5 parameters in results contain transactors
    filteredresults = filter(x -> any(occursin(t, x.TREATMENTS) for t in trans_actors), results)
    @test size(parameters, 1) == p*size(filteredresults, 1) == 25
    expected_variants = Dict(
        "RSID_103" => ["RSID_110", "RSID_10", "RSID_79", "RSID_6", "RSID_179"],
        "RSID_198" => ["RSID_20", "RSID_58", "RSID_197", "RSID_120", "RSID_65"]
    )
    for (initial_param_index, row) in enumerate(eachrow(filteredresults))
        bqtl, eqtl = NegativeControl.split_string(row.TREATMENTS)
        for new_param_index in 1:p
            Ψ = parameters[5*(initial_param_index - 1) + new_param_index]
            @test Ψ.target == Symbol(row.TARGET)
            @test Ψ.confounders == Symbol.(NegativeControl.split_string(row.CONFOUNDERS))
            @test Ψ.covariates == NegativeControl.getcovariates(row.COVARIATES)
            @test Ψ isa IATE
            @test keys(Ψ.treatment)[1] == Symbol(bqtl)
            @test keys(Ψ.treatment)[2] == Symbol(expected_variants[eqtl][new_param_index])
        end
    end
    
end

@testset "Test generate_random_variants_parameters_and_dataset" begin
    outdir = "resultdir"
    ispath(outdir) || mkdir(outdir)
    parsed_args = Dict(
        "p" => 5,
        "results" => joinpath("data", "summary.csv"),
        "trans-actors-prefix" => joinpath("data", "trans_act"),
        "bgen-prefix" => joinpath("data", "bgen", "ukb"),
        "outdir" => outdir,
        "pval-col" => "PVALUE",
        "pval-threshold" => 0.05,
        "verbosity" => 0,
        "reltol" => 0.05,
        "rng" => 123,
        "chunksize" => 15
    )  
    generate_random_variants_parameters_and_dataset(parsed_args)
    # Clean
    params_1 = deserialize(joinpath(outdir, "random_variants_param_1.bin"))
    params_2 = deserialize(joinpath(outdir, "random_variants_param_2.bin"))
    parameters = vcat(params_1, params_2)
    @test size(parameters, 1) == 25
    expected_mapped_variants = Dict(
        :RSID_110 => 0, 
        :RSID_10 => 0,
        :RSID_79 => 0,
        :RSID_179 => 0,
        :RSID_142 => 0,
        :RSID_20 => 0,
        :RSID_58 => 0,
        :RSID_120 => 0,
        :RSID_65 => 0,
        :RSID_127 => 0,
    )
    expected_targets = Dict(
        Symbol("K56 Paralytic ileus and intestinal obstruction without hernia") => 0,
        Symbol("E20-E35 Disorders of other endocrine glands") => 0,
        Symbol("O68 Labour and delivery complicated by foetal stress [distress]") => 0,
        Symbol("High light scatter reticulocyte percentage") => 0
    )
    for Ψ in parameters
        expected_mapped_variants[keys(Ψ.treatment)[2]] += 1
        @test Ψ.confounders == [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
        expected_targets[Ψ.target] += 1
    end
    @test expected_targets == Dict(
        Symbol("K56 Paralytic ileus and intestinal obstruction without hernia") => 5,
        Symbol("E20-E35 Disorders of other endocrine glands") => 5,
        Symbol("O68 Labour and delivery complicated by foetal stress [distress]") => 5,
        Symbol("High light scatter reticulocyte percentage") => 10
    )
    @test expected_mapped_variants == Dict(
        :RSID_179 => 3,
        :RSID_10  => 3,
        :RSID_58  => 2,
        :RSID_120 => 2,
        :RSID_65  => 2,
        :RSID_127 => 2,
        :RSID_79  => 3,
        :RSID_142 => 3,
        :RSID_110 => 3,
        :RSID_20  => 2
    )
    rm(outdir, force=true, recursive=true)
end

end

true