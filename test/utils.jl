module TestUtils

using Test
using NegativeControl

@testset "Test misc utils" begin
    # Test split_string
    @test NegativeControl.split_string("AC_&_CC") == ["AC", "CC"]

    # Test get confounders
    @test NegativeControl.getconfounders("PC1_&_PC2_&_PC3_&_PC4_&_Age") == [:Age]
    @test NegativeControl.getconfounders("PC1_&_PC2_&_PC3_&_PC5") == []
end

@testset "Test retrieve_significant_results" begin
    results_file = joinpath("data", "summary.csv")
    # The last is an ATE and filtered
    results = NegativeControl.retrieve_significant_results(results_file; threshold=:PVALUE => 1)
    @test size(results) == (8, 15)
    # With another P-value column constraint
    results = NegativeControl.retrieve_significant_results(results_file; threshold=:ADJUSTED_PVALUE => 0.9)
    @test size(results) == (7, 15)
end

end
