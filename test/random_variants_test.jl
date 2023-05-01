module RandomVariantsTests

using Test
using NegativeControl
using BGEN
using Serialization
using StableRNGs
using HTTP

@testset "Test is_in_regulatory_region" begin
    v = deserialize(joinpath("data", "variants", "not_annotated_variant.bin"))
    @test NegativeControl.is_in_regulatory_region(v) == false

    v = deserialize(joinpath("data", "variants", "regulatory_variant.bin"))
    @test NegativeControl.is_in_regulatory_region(v) == true

    v = deserialize(joinpath("data", "variants", "failing_variant.bin"))
    @test_throws HTTP.Exceptions.StatusError NegativeControl.is_in_regulatory_region(v)
end

@testset "Test find_maf_matching_random_variants" begin
    b = Bgen(
        BGEN.datadir("example.8bits.bgen"); 
        sample_path=BGEN.datadir("example.sample"), 
        idx_path=BGEN.datadir("example.8bits.bgen.bgi")
        )
    target_rsid = "RSID_103"
    variants = NegativeControl.find_maf_matching_random_variants!(
        b, target_rsid; 
        p=10, rng=StableRNG(123), reltol=0.05
    )
    @test size(variants, 1) == 10
    for v in variants
        @test rsid(v) != target_rsid
    end
end


end

true