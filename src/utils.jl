
"""
    instantiate_dataset(path::String)

Returns a DataFrame wrapper around a dataset, either in CSV format.
"""
instantiate_dataset(path::String) =
    endswith(path, ".csv") ? CSV.read(path, DataFrame, ntasks=1) : DataFrame(Arrow.Table(path))

BGEN.rsid(s::Symbol) = s

treatment_variables(Ψ::ComposedEstimand) =
    unique(vcat((treatment_variables(arg) for arg ∈ Ψ.args)...))

treatment_variables(Ψ) = collect(keys(Ψ.treatment_values))

outcome_variables(Ψ) = [Ψ.outcome]

outcome_variables(Ψ::ComposedEstimand) = 
    unique(vcat((outcome_variables(arg) for arg ∈ Ψ.args)...))

"""
    getconfounders(v)

Split string and remove principal components from the list.
"""
getconfounders(v) = Symbol.(filter(x -> !occursin(r"^PC[0-9]*$", x), split_string(v)))

default_statistical_test(Ψ̂; threshold=0.05) = pvalue(OneSampleTTest(Ψ̂)) < threshold

default_statistical_test(Ψ̂::TMLE.ComposedEstimate; threshold=0.05) = 
    length(Ψ̂.estimate) > 1 ? pvalue(TMLE.OneSampleHotellingT2Test(Ψ̂)) < threshold : pvalue(TMLE.OneSampleTTest(Ψ̂)) < threshold

is_significant(Ψ̂; threshold=0.05) = 
    default_statistical_test(Ψ̂; threshold=threshold)

function read_significant_from_hdf5(filename; threshold=0.05)
    jldopen(filename) do io
        return mapreduce(vcat, keys(io)) do key
            [Ψ̂.estimand for Ψ̂ ∈ io[key] if is_significant(Ψ̂; threshold=threshold)]
        end
    end
end

function read_significant_from_jls(filename; threshold=0.05)
    results = []
    open(filename) do io
        while !eof(io)
            Ψ̂ = deserialize(io)
            if is_significant(Ψ̂, threshold=threshold)
                push!(results, Ψ̂.estimand)
            end
        end
    end
    return results
end

read_significant_from_json(filename; threshold=0.05) = 
    [Ψ̂.estimand for Ψ̂ ∈ TMLE.read_json(filename) if is_significant(Ψ̂; threshold=threshold)]

function read_significant_results(filename; threshold=0.05)
    results = if endswith(filename, "hdf5")
        read_significant_from_hdf5(filename; threshold=threshold)
    elseif endswith(filename, "jls")
        read_significant_from_jls(filename; threshold=threshold)
    elseif endswith(filename, "json")
        read_significant_from_json(filename; threshold=threshold)
    else
        throw(ArgumentError("Unupported estimate file format: $filepath"))
    end
    return results
end

"""
    group_parameters(parameters; min_size=1)

Tries to group parameters together in optimal ordering 
by group of size approximately greater than min_size.
"""
function write_parameter_files(outdir, parameters, chunksize; prefix="permutation_estimands_")
    for (index, batch) in enumerate(Iterators.partition(parameters, chunksize))
        serialize(joinpath(outdir, string(prefix, index, ".yaml")), Configuration(estimands=batch))
    end
end

function unique_treatments(estimands)
    treatments = Set{String}()
    for Ψ in estimands
        union!(treatments, string.(treatment_variables(Ψ)))
    end
    return treatments
end