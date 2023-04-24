#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.RESULTS_FILE = ""
params.RANDOM_VARIANTS = ""

process GeneratePermutationTestsData {
    container "olivierlabayle/neg-cointrol:main"
    publishDir "$params.OUTDIR/parameters", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR/tmle_inputs", mode: 'symlink', pattern: "*.csv"
    label "bigmem"
    
    input:
        results_file

    output:
        path "dataset.csv", emit: dataset
        path "final.*.yaml", emit: parameters

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        """
}

workflow {
    results_file = Channel.value(file(params.RESULTS_FILE))

    dataset, param_files = GeneratePermutationTestsData(results_file)
    RunPermutationTests(dataset, param_files)


}