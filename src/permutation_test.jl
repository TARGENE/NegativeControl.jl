
permuted_name(colname) = string(colname, "_permuted")
split_string(s) = split(s, "_&_")

"""
Retrieve significant results defined by a threshold `Pair` `colname => threshold` 
from a set of estimation results given by `filepath` 
"""
function retrieve_significant_results(filepath; threshold=:PVALUE => 0.05)
    data = CSV.read(filepath, DataFrame)
    pvalcol, pval = threshold
    return filter(x -> x[pvalcol] < pval && x.PARAMETER_TYPE =="IATE" , data)
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
function build_negative_control_dataset(datafilepath, results; rng=StableRNG(123))
    data = CSV.read(datafilepath, DataFrame)
    treatment_cols = unique(reduce(vcat, [split(x, "_&_") for x in results[!, :TREATMENTS]]))
    target_cols = unique(results.TARGET)
    for col in vcat(treatment_cols, target_cols)
        make_permuted_col!(data, col; rng=rng)
    end
    return data
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
        covariates  = Symbol.(split_string(param_row.COVARIATES))
    )
end

"""
New parameters are generated by permuting treatment and target columns
at each order <= interaction order.
"""
function make_permutation_parameters(results)
    parameters = IATE[]
    for (treatments, treatment_grp) in pairs(groupby(results, :TREATMENTS))
        treatments = split_string.(treatments.TREATMENTS)   
        for (target, target_grp) in pairs(groupby(treatment_grp, :TARGET))
            target = target.TARGET
            for comb in combinations(vcat(treatments, target))
                new_treatment, new_target = permutation_setting(comb, treatments, target)
                parameters = vcat(
                    parameters, 
                    make_parameters(target_grp, new_target, new_treatment),
                )
            end
        end
    end
    return sort(parameters, by= x -> tuple(keys(x.treatment)..., x.target))
end

function run_permutation_test(parsed_args)
    resultdir = parsed_args["resultdir"]
    outdir = parsed_args["outdir"]
    pval_col = parsed_args["pval-col"]
    pval_threshold = parsed_args["pval-threshold"]
    estimator_file = parsed_args["estimator-file"]
    verbosity = parsed_args["verbosity"]
    limit = parsed_args["limit"]
    rng_int = parsed_args["rng"]

    ispath(outdir) || mkpath(outdir)
    verbosity > 0 && @info "Loading TarGene run summary.csv"
    results = NegativeControl.retrieve_significant_results(joinpath(resultdir, "summary.csv"), threshold=pval_col => pval_threshold)
    verbosity > 0 && @info string(size(results, 1), " parameters satisfying the threshold.")
    verbosity > 0 && @info "Generating permutation parameters."
    parameters = NegativeControl.make_permutation_parameters(results)
    if limit !== nothing
        parameters = parameters[1:limit]
    end
    verbosity > 0 && @info string(size(parameters, 1), " parameters will be estimated.")
    verbosity > 0 && @info "Building permutation dataset"
    dataset = NegativeControl.build_negative_control_dataset(joinpath(resultdir, "tmle_inputs", "final.data.csv"), results; rng=StableRNG(rng_int))

    treatment_cols = unique(Iterators.flatten(keys(p.treatment) for p in parameters))
    covariate_cols = unique(Iterators.flatten(p.covariates for p in parameters))
    confounder_cols = unique(Iterators.flatten(p.confounders for p in parameters))
    TargetedEstimation.make_categorical!(dataset, treatment_cols, infer_ordered=true)
    
    # Confounders and Covariates are converted to Float64
    TargetedEstimation.make_float!(dataset, vcat(covariate_cols, confounder_cols))
    
    # Retrieve TMLE specifications
    tmle_spec = TargetedEstimation.tmle_spec_from_yaml(estimator_file)
    csv_io = TargetedEstimation.initialize_csv_io(joinpath(outdir, "summary"))
    cache = TMLECache(dataset)
    for Ψ in parameters
        target = Ψ.target
        targetisbinary = TargetedEstimation.isbinarytarget(dataset[!, target])
        targetisbinary && TargetedEstimation.make_categorical!(dataset, target)
        η_spec = TargetedEstimation.nuisance_spec_from_target(tmle_spec, targetisbinary, tmle_spec.cache)
        tmle_result, log = TargetedEstimation.try_tmle!(cache, Ψ, η_spec; verbosity=verbosity, threshold=tmle_spec.threshold)
        # Append CSV result for target
        TargetedEstimation.append_csv(csv_io, (PARAMETER=[Ψ],), [tmle_result], [log])
    end

    verbosity > 0 && @info "Done."
end
