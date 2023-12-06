module NegativeControl

using CSV 
using DataFrames
using StableRNGs
using Random
using Combinatorics
using TMLE
using Arrow
using Serialization
using BGEN
using JSON
using HTTP
using Statistics
using JLD2

include("utils.jl")
include("permutation_test.jl")
include("random_variants_test.jl")


export generate_permutation_parameters_and_dataset, generate_random_variants_parameters_and_dataset

end
