using ArgParse
using NegativeControl

function parse_commandline()
    s = ArgParseSettings(
        description = "Random variant parameters generation",
        commands_are_required = false,
        version = "0.1",
        add_version = true)

    @add_arg_table s begin
        "trans-actors"
            help = "Trans acting factors in csv format"
            required = true
        "results"
            help = "CSV file containing TarGene results."
            required = true
        "bgen-prefix"
            help = "BGEN chromosome files prefix"
            required = true
        "--outdir"
            help = "Where the data will be generated."
            arg_type = String
            default = pwd()
        "--p"
            help = "Number of random variants per trans-actor"
            arg_type = Int
            default = 10
        "--pval-col"
            help = "The p-value column to search for significant results"
            default = "TMLE_PVALUE"
            arg_type = String
        "--pval-threshold"
            help = "The p-value threshold for significant results calling"
            default = 0.05
            arg_type = Float64
        "--reltol"
            help = "Relative tolerance between a trans-actor MAF and its mapped variants"
            default = 0.05
            arg_type = Float64
        "--rng"
            help = "The random seed for permutations"
            default = 123
            arg_type = Int
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

generate_random_variants_parameters_and_dataset(parsed_args)


