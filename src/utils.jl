
BGEN.rsid(s::AbstractString) = s

split_string(s) = split(s, "_&_")

getcovariates(v::Missing) = Symbol[]
getcovariates(v) = Symbol.(split_string(v))

"""
Retrieve significant results defined by a threshold `Pair` `colname => threshold` 
from a set of estimation results given by `filepath` 
"""
function retrieve_significant_results(filepath; threshold=:PVALUE => 0.05)
    data = CSV.read(filepath, DataFrame)
    pvalcol, pval = threshold
    return filter(x -> x[pvalcol] !== missing && x[pvalcol] < pval && x.PARAMETER_TYPE =="IATE" , data)
end

"""
    group_parameters(parameters; min_size=1)

Tries to group parameters together in optimal ordering 
by group of size approximately greater than min_size.
"""
function write_parameter_files(outdir, parameters, chunksize; prefix="permutation_param_")
    for (index, param_group) in enumerate(Iterators.partition(parameters, chunksize))
        serialize(joinpath(outdir, string(prefix, index, ".bin")), param_group)
    end
end

function unique_treatments(results::DataFrame)
    treatments = Set{String}()
    for treatment_string in results.TREATMENTS
        push!(treatments, split_string(treatment_string)...)
    end
    return treatments
end