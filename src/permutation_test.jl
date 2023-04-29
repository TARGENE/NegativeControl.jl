
permuted_name(colname) = string(colname, "_permuted")
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

function permutation_setting(comb, treatments, target)
    new_treatments = []
    for treatment in treatments
        if treatment ∈ comb
            push!(new_treatments, permuted_name(treatment))
        else
            push!(new_treatments, treatment)
        end
    end
    new_target = target ∈ comb ? permuted_name(target) : target
    return new_treatments, new_target
end

make_permuted_col!(data, col; rng=StableRNG(123)) = 
    data[!, permuted_name(col)] = shuffle(rng, data[!, col])

"""
Rewrite the TMLE source datafile with additional permuted treatment columns
"""
function write_negative_control_dataset(outfile, datafilepath, results; rng=StableRNG(123))
    data = CSV.read(datafilepath, DataFrame)
    treatment_cols = unique(reduce(vcat, [split(x, "_&_") for x in results[!, :TREATMENTS]]))
    target_cols = unique(results.TARGET)
    for col in vcat(treatment_cols, target_cols)
        make_permuted_col!(data, col; rng=rng)
    end
    Arrow.write(outfile, data)
end

function parse_case_control(c)
    c_string = split_string(c)
    final_c = Vector(undef, size(c_string, 1))
    for index in eachindex(final_c)
        c_ = tryparse(Int, c_string[index])
        final_c[index] = c_ === nothing ? c_string[index] : c_
    end
    return final_c
end

function make_parameters(df, target, treatments)
    params = Vector{IATE}(undef, size(df, 1))
    for (index, param_row) in enumerate(eachrow(df))
        params[index] = make_parameter(param_row, target, treatments)
    end
    return params
end

function make_parameter(param_row, target, treatments)
    case = parse_case_control(param_row.CASE)
    control = parse_case_control(param_row.CONTROL)
    treatment = NamedTuple{Tuple(Symbol.(treatments))}([(case=cs, control=ct) for (cs, ct) in zip(case, control)])
    return IATE(
        target      = Symbol(target),
        treatment   = treatment,
        confounders = Symbol.(split_string(param_row.CONFOUNDERS)),
        covariates  = covariates(param_row.COVARIATES) 
    )
end

"""
New parameters are generated by permuting treatment and target columns
at each order <= interaction order.
"""
function make_permutation_parameters(results, orders)
    parameters = IATE[]
    for (treatments, treatment_grp) in pairs(groupby(results, :TREATMENTS))
        treatments = split_string.(treatments.TREATMENTS)   
        for (target, target_grp) in pairs(groupby(treatment_grp, :TARGET))
            target = target.TARGET
            for order in orders
                combs = combinations(vcat(treatments, target), order) 
                for comb in combs
                    new_treatment, new_target = permutation_setting(comb, treatments, target)
                    parameters = vcat(
                        parameters, 
                        make_parameters(target_grp, new_target, new_treatment),
                    )
                end
            end
        end
    end
    return optimize_ordering(parameters)
end

"""
    group_parameters(parameters; min_size=1)

Tries to group parameters together in optimal ordering 
by group of size approximately greater than min_size.
"""
function write_parameter_files(outdir, parameters, chunksize)
    for (index, param_group) in enumerate(Iterators.partition(parameters, chunksize))
        serialize(joinpath(outdir, string("param_", index, ".bin")), param_group)
    end
end

function generate_permutation_parameters_and_dataset(parsed_args)
    # Parsing Arguments
    datafile = parsed_args["dataset"]
    resultsfile = parsed_args["results"]
    outdir = parsed_args["outdir"]
    pval_col = parsed_args["pval-col"]
    pval_threshold = parsed_args["pval-threshold"]    
    verbosity = parsed_args["verbosity"]
    orders = parse.(Int, split(parsed_args["orders"], ","))
    limit = parsed_args["limit"]
    rng_int = parsed_args["rng"]
    chunksize = parsed_args["chunksize"]

    # Generating Permutation Parameters
    verbosity > 0 && @info string("Retrieving significant parameters.")
    results = retrieve_significant_results(resultsfile, threshold=pval_col => pval_threshold)
    verbosity > 0 && @info string(size(results, 1), " parameters satisfying the threshold.")
    verbosity > 0 && @info "Generating permutation parameters."
    parameters = make_permutation_parameters(results, orders)
    if limit !== nothing
        parameters = parameters[1:limit]
    end
    verbosity > 0 && @info string(size(parameters, 1), " parameters will be estimated.")
    write_parameter_files(outdir, parameters, chunksize)
    
    # Building Permutation Dataset
    verbosity > 0 && @info "Building permutation dataset"
    write_negative_control_dataset(
        joinpath(outdir, "permutation_dataset.arrow"), 
        datafile, 
        results; 
        rng=StableRNG(rng_int)
    )
    verbosity > 0 && @info "Done."
end
