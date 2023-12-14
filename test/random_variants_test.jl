module RandomVariantsTests

using Test
using NegativeControl
using BGEN
using Serialization
using StableRNGs
using HTTP
using Statistics
using TMLE

TESTDIR = joinpath(pkgdir(NegativeControl), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test same_maf and isSNP" begin
    b = Bgen(
        BGEN.datadir("example.8bits.bgen"); 
        sample_path=BGEN.datadir("example.sample"), 
        idx_path=BGEN.datadir("example.8bits.bgen.bgi")
    )
    v = variant_by_rsid(b, "RSID_198")
    @test mean(minor_allele_dosage!(b, v)) ≈ 0.48411763 atol=1e-5
    @test NegativeControl.same_maf(b, v, 0.51; reltol=0.05) == false
    @test NegativeControl.same_maf(b, v, 0.50; reltol=0.05) == true

    @test NegativeControl.isSNP(v)
    v.alleles = ["AG", "GG"]
    @test !NegativeControl.isSNP(v)
end

@testset "Test get_new_treatment_setting" begin
    b = Bgen(
        BGEN.datadir("example.8bits.bgen"); 
        sample_path=BGEN.datadir("example.sample"), 
        idx_path=BGEN.datadir("example.8bits.bgen.bgi")
    )
    v₁ = variant_by_rsid(b, "RSID_198")
    minor_allele_dosage!(b, v₁)
    v₂ = variant_by_rsid(b, "RSID_199")
    minor_allele_dosage!(b, v₂)

    @test minor_allele(v₁) == "A"
    @test major_allele(v₁) == "G"
    @test minor_allele(v₂) == "G"
    @test major_allele(v₂) == "A"
    origin_setting = (control="AA", case="AG")
    @test NegativeControl.get_new_treatment_setting(v₁, v₂, origin_setting) == (control="GG", case="GA")
    origin_setting = (control="AA", case="GG")
    @test NegativeControl.get_new_treatment_setting(v₁, v₂, origin_setting) == (control="GG", case="AA")
    origin_setting = (control="AG", case="AA")
    @test NegativeControl.get_new_treatment_setting(v₁, v₂, origin_setting) == (control="GA", case="GG")
end

@testset "Test is_in_regulatory_region" begin
    v = deserialize(joinpath(TESTDIR, "data", "variants", "not_annotated_variant.bin"))
    @test NegativeControl.is_in_regulatory_region(v) == false

    v = deserialize(joinpath(TESTDIR, "data", "variants", "regulatory_variant.bin"))
    @test NegativeControl.is_in_regulatory_region(v) == true

    v = deserialize(joinpath(TESTDIR, "data", "variants", "failing_variant.bin"))
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
        p=10, rng=StableRNG(123), reltol=0.05, verbosity=0
    )
    @test size(variants, 1) == 10
    for v in variants
        @test rsid(v) != "RSID_103"
    end

    # All transactors
    p = 10
    bgen_prefix = joinpath(TESTDIR, "data", "bgen", "ukb")
    trans_actors = Set(["RSID_103", "RSID_198"])
    reltol = 0.1
    variant_map = NegativeControl.find_maf_matching_random_variants(
        trans_actors, 
        bgen_prefix,
        trans_actors; 
        p=p, rng=StableRNG(123), reltol=reltol, verbosity=0
    )

    # Check criteria
    expected_mapped_rsids = Dict(
        :RSID_103 => [
            "RSID_110", "RSID_10", "RSID_79", "RSID_6", "RSID_179", 
            "RSID_142", "RSID_42", "RSID_106", "RSID_61", "RSID_57"
            ], 
        :RSID_198 => [
            "RSID_20", "RSID_58", "RSID_197", "RSID_120", "RSID_65", 
            "RSID_127", "RSID_98", "RSID_188", "RSID_165", "RSID_158",
            ]
    )
    for (key, mapped_rsids) in expected_mapped_rsids
        variant = variant_map[key][1]
        mapped_variants = variant_map[key][2]
        @test Symbol(rsid(variant)) == key
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
        p=p, rng=StableRNG(123), reltol=reltol, verbosity=0
    )
end

@testset "Test make_random_variants_parameters" begin
    p = 5
    reltol = 0.1
    rng = StableRNG(123)
    bgen_prefix = joinpath(TESTDIR, "data", "bgen", "ukb")
    trans_actors = Set(["RSID_103", "RSID_104"])
    
    estimands = [nt.TMLE.estimand for nt ∈ make_estimates()[1:end-1]]
    variant_map = NegativeControl.find_maf_matching_random_variants(
        trans_actors, 
        bgen_prefix,
        trans_actors; 
        p=p, rng=rng, reltol=reltol, verbosity=0
    )
    new_estimands = NegativeControl.make_random_variants_estimands(estimands, variant_map, p=p)
    # 5 new estimands for each of the 4 input estimands
    @test length(new_estimands) == 20
    # Check 5 new estimands generated from estimand 1
    from_1 = new_estimands[1:5]
    @test all(Ψ isa TMLE.StatisticalIATE for Ψ ∈ from_1)
    expected_new_rsids = [Symbol(rsid(v)) for v ∈ variant_map[:RSID_103][2]]
    new_rsids = [keys(Ψ.treatment_values)[1]  for Ψ ∈ from_1]
    @test expected_new_rsids == new_rsids
    # Check 5 new estimands generated from estimand 2
    from_2 = new_estimands[6:10]
    @test all(Ψ isa TMLE.StatisticalIATE for Ψ ∈ from_2)
    new_rsids = [keys(Ψ.treatment_values)[1]  for Ψ ∈ from_2]
    @test expected_new_rsids == new_rsids
    # Check 5 new estimands generated from estimand 3 (composed estimand)
    for Ψ ∈ new_estimands[11:15]
        @test Ψ isa ComposedEstimand
        @test Ψ != estimands[3]
    end
    # Check 5 new estimands generated from estimand 4
    from_4 = new_estimands[16:20]
    @test all(Ψ isa TMLE.StatisticalATE for Ψ ∈ from_4)
    expected_new_rsids = [Symbol(rsid(v)) for v ∈ variant_map[:RSID_104][2]]
    new_rsids = [keys(Ψ.treatment_values)[1] for Ψ ∈ from_4]
    @test expected_new_rsids == new_rsids
end

@testset "Test generate_random_variants_parameters_and_dataset" begin
    estimates = make_estimates()
    save(estimates)
    parsed_args = Dict(
        "p" => 5,
        "results" => "tmle_output.hdf5",
        "trans-actors-prefix" => joinpath(TESTDIR, "data", "trans_act"),
        "bgen-prefix" => joinpath(TESTDIR, "data", "bgen", "ukb"),
        "out" => "random_variants_parameters.jls",
        "pval-threshold" => 0.05,
        "estimator-key" => "TMLE",
        "verbosity" => 0,
        "reltol" => 0.05,
        "rng" => 123,
    )  
    generate_random_variants_parameters_and_dataset(parsed_args)
    
    # Only RSID_103 is a trans actor, leading to 3*5=15 estimands
    new_estimands = deserialize(parsed_args["out"]).estimands
    @test length(new_estimands) == 15
    # Check counts
    counts = Dict(
        "Composed" => 0,
        "IATE" => 0,
    )
    for Ψ ∈ new_estimands
        if Ψ isa ComposedEstimand
            counts["Composed"] += 1
        elseif Ψ isa TMLE.StatisticalIATE
            counts["IATE"] += 1
        end
    end
    @test counts == Dict(
        "Composed" => 5,
        "IATE" => 10
    )

    # Clean
    clean()
    rm(parsed_args["out"])
end

end

true