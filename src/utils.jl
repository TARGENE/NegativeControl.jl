
"""
    instantiate_dataset(path::String)

Returns a DataFrame wrapper around a dataset, either in CSV format.
"""
instantiate_dataset(path::String) =
    endswith(path, ".csv") ? CSV.read(path, DataFrame, ntasks=1) : DataFrame(Arrow.Table(path))

    
BGEN.rsid(s::AbstractString) = s

split_string(s) = split(s, "_&_")

getcovariates(v::Missing) = Symbol[]
getcovariates(v) = Symbol.(split_string(v))

"""
    getconfounders(v)

Split string and remove principal components from the list.
"""
getconfounders(v) = Symbol.(filter(x -> !occursin(r"^PC[0-9]*$", x), split_string(v)))


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
        parameters_to_yaml(joinpath(outdir, string(prefix, index, ".yaml")), param_group)
    end
end

function unique_treatments(results::DataFrame)
    treatments = Set{String}()
    for treatment_string in results.TREATMENTS
        push!(treatments, split_string(treatment_string)...)
    end
    return treatments
end