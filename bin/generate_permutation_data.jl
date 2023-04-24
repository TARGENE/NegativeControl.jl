using ArgParse
using NegativeControl

function parse_commandline()
    s = ArgParseSettings(
        description = "Permutation tests data generation",
        commands_are_required = false,
        version = "0.1",
        add_version = true)

    @add_arg_table s begin
        "dataset"
            help = "Path to tmle input dataset file (.csv|.arrow)"
            required = true
        "results"
            help = "CSV file containing TarGene results."
            required = true
        "--outdir"
            help = "Where the data will be generated."
            arg_type = String
            default = pwd()
        "--pval-col"
            help = "The p-value column to search for significant results"
            default = "TMLE_PVALUE"
            arg_type = String
        "--pval-threshold"
            help = "The p-value threshold for significant results calling"
            default = 0.05
            arg_type = Float64
        "limit"
            help = "The max number of permutation parameters to be generated"
            default = nothing
            arg_type = Int
        "rng"
            help = "The random seed for permutations"
            default = 123
            arg_type = Int
        "orders"
            help = "A comma separated set of combination orders e.g. 1,2,3"
            default = "1"
            arg_type = String
        "--chunksize"
            help = "Results will be appended to outfiles every chunk"
            default = 100
            arg_type = Int
        "--verbosity", "-v"
            help = "Verbosity level"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

generate_permutation_parameters_and_dataset(parsed_args)


