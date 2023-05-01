
split_string(s) = split(s, "_&_")

covariates(v::Missing) = Symbol[]
covariates(v) = Symbol.(split_string(v))

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
function write_parameter_files(outdir, parameters, chunksize)
    for (index, param_group) in enumerate(Iterators.partition(parameters, chunksize))
        serialize(joinpath(outdir, string("permutation_param_", index, ".bin")), param_group)
    end
end