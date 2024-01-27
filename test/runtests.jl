using Test
using NegativeControl

TESTDIR = joinpath(pkgdir(NegativeControl), "test")

@test include(joinpath(TESTDIR, "utils.jl"))
@test include(joinpath(TESTDIR, "permutation_test.jl"))
@test include(joinpath(TESTDIR, "random_variants_test.jl"))