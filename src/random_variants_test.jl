const ENSEMBL_URL = "https://rest.ensembl.org"
const REGULATORY_ELEMENTS = Set(["exon"])
const CHR_REG = r"chr[1-9]+"

NotEnoughMatchingVariantsError(rsid, p, reltol) = 
    ArgumentError(string("Could not find ", p, 
    " random matching variants for ", rsid, ", consider decreasing the MAF reltol (current: ", reltol, ")"
    ))


read_snps_from_csv(path::String) = unique(CSV.read(path, DataFrame; select=[:ID, :CHR]), :ID)

function trans_actors_from_prefix(trans_actors_prefix::AbstractString)
    dir, prefix = splitdir(trans_actors_prefix)
    dir_ = dir == "" ? "." : dir
    trans_actors = DataFrame[]
    for filename in readdir(dir_)
        if startswith(filename, prefix)
            push!(trans_actors, read_snps_from_csv(joinpath(dir, filename)))
        end
    end
    return Set(vcat((df.ID for df in trans_actors)...))
end

"""
    control_case(variant::Variant, origin_variant::Variant, origin_case::Int, origin_control::Int)

If the allele representation is an Integer then it is kep the same.
"""
function control_case(variant::Variant, origin_variant::Variant, origin_case::Int, origin_control::Int)
    return (control=origin_control, case=origin_case)
end

"""
    control_case(variant::Variant, origin_variant::Variant, origin_case::AbstractString, origin_control::AbstractString)

If the allele representation is a String then new corresponding alleles are inferred for the mapped variant.
"""
function control_case(variant::Variant, origin_variant::Variant, origin_case::AbstractString, origin_control::AbstractString)
    allele_map = Dict(
        minor_allele(origin_variant) => minor_allele(variant),
        major_allele(origin_variant) => major_allele(variant)
    )
    control = join(allele_map[string(a)] for a in origin_control)
    case = join(allele_map[string(a)] for a in origin_case)
    return (control=control, case=case)
end

"""
    update_treatment_setting!(treatment_setting, variant::Variant, origin_variant::Variant, origin_case, origin_control)

When the origin_variant is a `Variant` it means it is a "mapped" transactor and
new alleles must be inferred.
"""
function update_treatment_setting!(treatment_setting, variant::Variant, origin_variant::Variant, origin_case, origin_control)
    push!(
        treatment_setting, 
        control_case(variant, origin_variant, origin_case, origin_control)
        )
end

"""
    update_treatment_setting!(treatment_setting, variant::Variant, origin_variant::AbstractString, origin_case, origin_control)

When the origin_variant is an `AbstractString` it means it is a "non-mapped" treatment and
treatment settings are kept identical.
"""
function update_treatment_setting!(treatment_setting, variant::AbstractString, origin_variant::AbstractString, origin_case, origin_control)
    push!(
        treatment_setting, 
        (control=origin_control, case=origin_case)
        )
end

function is_numbered_chromosome_file(filename, prefix)
    if occursin(prefix, filename) && endswith(filename, "bgen")
        regexp_match = match(CHR_REG, filename)
        if regexp_match !== nothing
            return true
        end
    end
    return false
end

function read_bgen(filepath)
    sample_filepath = string(filepath[1:end-4], "sample")
    idx_filepath = string(filepath, ".bgi")
    return Bgen(filepath, sample_path=sample_filepath, idx_path=idx_filepath)
end

"""
    overlaps(v::Variant)

Checks whether a variant is in a regulatory region as defined by ENSEMBL. nm,./\
"""
function is_in_regulatory_region(v::Variant)
    chr = parse(Int, chrom(v))
    loc = pos(v)
    ext = string("/overlap/region/human/", chr, ":", loc, "-", loc, "?feature=regulatory;feature=exon;feature=other_regulatory")
    url = string(ENSEMBL_URL, ext)
    resp = HTTP.get(url, headers=Dict( "Content-Type" => "application/json"))
    if length(JSON.parse(String(resp.body))) == 0
        return false
    else
        return true
    end
end

function same_maf(b::Bgen, variant::Variant, maf::AbstractFloat; reltol=0.05)
    variant_maf = mean(minor_allele_dosage!(b, variant))
    return abs(maf - variant_maf)/maf <= reltol
end

function matches_constraints!(b::Bgen, candidate_variant::Variant, target_maf::AbstractFloat, rs_ids::Set{<:AbstractString}; reltol=0.05)
    rsid(candidate_variant) ∈ rs_ids && return false
    if same_maf(b, candidate_variant, target_maf; reltol=reltol)
        isreg = false
        try
            isreg = is_in_regulatory_region(candidate_variant)
        catch
            return false
        end
        return !isreg
    end
    return false
end

"""
find_maf_matching_random_variants!(b::Bgen, variant::Variant, rs_ids::Set{<:AbstractString}; p=10, rng=StableRNG(123), reltol=0.05)

For the given variant find `p` random variants in the BGEN file
that match the variant's MAF (up to the relative `tol`) and are not in:

- (i) exons
- (ii) promoters
- (iii) enhancers
- (iv) `rs_ids`

"""
function find_maf_matching_random_variants!(b::Bgen, variant::Variant, rs_ids::Set{<:AbstractString}; p=10, rng=StableRNG(123), reltol=0.05)
    target_maf = mean(minor_allele_dosage!(b, variant))
    all_rsids = shuffle(rng, rsids(b))
    matching_variants = Vector{Variant}(undef, p)
    index = 1
    for candidate_rs_id in all_rsids
        candidate_variant = variant_by_rsid(b, candidate_rs_id)
        if matches_constraints!(b, candidate_variant, target_maf, rs_ids; reltol=reltol)
            matching_variants[index] = candidate_variant
            index += 1
            index > p && return matching_variants
        end
    end
    throw(NotEnoughMatchingVariantsError(rsid(variant), p, reltol))
end

function find_maf_matching_random_variants(trans_actors::Set{<:AbstractString}, bgen_prefix::AbstractString, all_rsids::Set{<:AbstractString}; p=10, rng=StableRNG(123), reltol=0.05)
    chr_dir_, prefix_ = splitdir(bgen_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    variant_map = Dict()
    remaining_transactors = deepcopy(trans_actors)
    for filename in readdir(chr_dir)
        length(remaining_transactors) == 0 && return variant_map
        if is_numbered_chromosome_file(filename, prefix_)
            b = read_bgen(joinpath(chr_dir_, filename))
            bgen_rsids = Set(rsids(b))
            for rs_id in remaining_transactors
                if rs_id ∈ bgen_rsids
                    variant = variant_by_rsid(b, rs_id)
                    minor_allele_dosage!(b, variant)
                    variant_map[rs_id] = (
                        variant,
                        find_maf_matching_random_variants!(
                        b, variant, all_rsids;
                        p=p, rng=rng, reltol=reltol)
                    )
                    pop!(remaining_transactors, rs_id)
                end
            end
        end
    end
    return variant_map
end


function make_random_variants_parameters(results, variant_map)
    parameters = TMLE.Parameter[]
    transactors = keys(variant_map)
    
    for row in eachrow(results)
        # At least one trans-actor in the parameter treatments to be processed
        if any(occursin(rs_id, row.TREATMENTS) for rs_id in transactors)
            confounders = Symbol.(split_string(row.CONFOUNDERS))
            target = Symbol(row.TARGET)
            covariates = getcovariates(row.COVARIATES)
            paramtype = getfield(TMLE, Symbol(row.PARAMETER_TYPE))

            control = split_string(row.CONTROL)
            case = split_string(row.CASE)
            if tryparse(Int, first(control)) !== nothing
                control = parse.(Int, control)
                case = parse.(Int, case)
            end

            origin_treatments_ids = split_string(row.TREATMENTS)
            treatments = [t ∈ transactors ? variant_map[t][2] : [t] for t in origin_treatments_ids]
            for treatment_vars ∈ Iterators.product(treatments...)
                treatment_setting = []
                for (index, treatment_var) in enumerate(treatment_vars)
                    origin_id = origin_treatments_ids[index]
                    origin_variant = haskey(variant_map, origin_id) ? variant_map[origin_id][1] : treatment_var
                    update_treatment_setting!(
                        treatment_setting, 
                        treatment_var, 
                        origin_variant, 
                        case[index], 
                        control[index]
                    )
                end
                treatment = NamedTuple{Tuple(Symbol(rsid(v)) for v in treatment_vars)}(treatment_setting)
                push!(
                    parameters,
                    paramtype(
                        target=target,
                        treatment=treatment,
                        confounders=confounders,
                        covariates=covariates
                    )
                )
            end
        end
    end
    return parameters
end


function generate_random_variants_parameters_and_dataset(parsed_args)
    resultsfile = parsed_args["results"]
    p = parsed_args["p"]
    rng = StableRNG(parsed_args["rng"])
    reltol = parsed_args["reltol"]
    trans_actors = trans_actors_from_prefix(parsed_args["trans-actors-prefix"])
    bgen_prefix = parsed_args["bgen-prefix"]
    pval_col = parsed_args["pval-col"]
    pval_threshold = parsed_args["pval-threshold"]
    out = parsed_args["out"]
    chunksize = parsed_args["chunksize"]
    verbosity = parsed_args["verbosity"]
    
    verbosity > 0 && @info string("Retrieving significant parameters.")
    results = retrieve_significant_results(resultsfile, threshold=pval_col => pval_threshold)
    all_rsids = unique_treatments(results)

    verbosity > 0 && @info string("Looking for random MAF matching variants for each trans-actor.")
    variant_map = find_maf_matching_random_variants(trans_actors, bgen_prefix, all_rsids; p=p, rng=rng, reltol=reltol)
    verbosity > 0 && @info string("Building new parameters.")
    parameters = make_random_variants_parameters(results, variant_map)
    optimize_ordering!(parameters)
    parameters_to_yaml(out, parameters)
    verbosity > 0 && @info string("Done.")
end