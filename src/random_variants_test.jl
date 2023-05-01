const ENSEMBL_URL = "https://rest.ensembl.org"
const REGULATORY_ELEMENTS = Set(["exon"])
const CHR_REG = r"chr[1-9]+"

NotEnoughMatchingVariantsError(rsid, p, reltol) = 
    ArgumentError(string("Could not find ", p, 
    " random matching variants for ", rsid, ", consider decreasing the MAF reltol (current: ", reltol, ")"
    ))

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

function matches_constraints!(b, candidate_variant, target_maf; reltol=0.05)
    candidate_maf = mean(minor_allele_dosage!(b, candidate_variant))
    if abs(target_maf - candidate_maf)/target_maf <= reltol
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
    find_maf_matching_random_variants(rsid, b::Bgen; p=10, rng=StableRNG(123))

For the given variant identified by its `rsid``, find `p` random variants in the BGEN file
that match the variant's MAF (up to the relative `tol`) and are not in:

- (i) exons
- (ii) promoters
- (iii) enhancers

"""
function find_maf_matching_random_variants!(b::Bgen, rs_id::String; p=10, rng=StableRNG(123), reltol=0.05)
    variant = variant_by_rsid(b, rs_id)
    target_maf = mean(minor_allele_dosage!(b, variant))
    all_rsids = shuffle(rng, rsids(b))
    matching_variants = Vector{Variant}(undef, p)
    index = 1
    for candidate_rs_id in all_rsids
        candidate_rs_id == rsid(variant) && continue
        candidate_variant = variant_by_rsid(b, candidate_rs_id)
        if matches_constraints!(b, candidate_variant, target_maf)
            matching_variants[index] = candidate_variant
            index += 1
            index > p && return matching_variants
        end
    end
    throw(NotEnoughMatchingVariantsError(rsid, p, reltol))
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

function find_maf_matching_random_variants!(bgen_prefix::String, trans_actors::Set{String}; p=10, rng=StableRNG(123), reltol=0.05)
    chr_dir_, prefix_ = splitdir(bgen_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    variant_map = Dict(rsid_ => Variant[] for rsid_ in trans_actors)
    for filename in readdir(chr_dir)
        length(trans_actors) == 0 && return variant_map
        if is_numbered_chromosome_file(filename, prefix_)
            b = read_bgen(joinpath(chr_dir_, filename))
            bgen_rsids = Set(rsids(b))
            for rs_id in trans_actors.ID
                if rs_id ∈ bgen_rsids
                    variant_map[rs_id] = find_maf_matching_random_variants!(
                        b, rs_id; 
                        p=p, rng=rng, reltol=reltol
                    )
                    pop!(trans_actors, rs_id)
                end
            end
        end
    end
    return variant_map
end


function make_random_variants_parameters(results, variant_map)
    
    for row in eachrow(results)

    end
end


function generate_random_variants_parameters_and_dataset(parsed_args)
    resultsfile = parsed_args["results"]
    p = parsed_args["p"]
    rng = parsed_args["rng"]
    reltol = parsed_args["reltol"]
    trans_actors = Set(CSV.read(parsed_args["trans-actors"], DataFrame).ID)
    bgen_prefix = parsed_args["bgen-prefix"]
    pval_col = parsed_args["pval-col"]
    pval_threshold = parsed_args["pval-threshold"]

    results = retrieve_significant_results(resultsfile, threshold=pval_col => pval_threshold)
    variant_map = find_maf_matching_random_variants!(bgen_prefix, trans_actors; p=p, rng=rng, reltol=reltol)
    new_parameters = make_random_variants_parameters(results, variant_map)
end