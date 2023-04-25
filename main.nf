#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Global Parameters

params.RESULTS_FILE = ""
params.RANDOM_VARIANTS = ""
params.VERBOSITY = 1

// Permutation Tests Parameters
params.MAX_PERMUTATION_TESTS = null
params.PVAL_COL = "ADJUSTED_PVALUE"
params.PVAL_THRESHOLD = 0.05
params.PERMUTATION_CHUNKSIZE
params.PERMUTATION_ORDERS = "1"
params.PERMUTATION_RNG = 123

process GeneratePermutationTestsData {
    container "olivierlabayle/negative-cointrols:initial_pipeline"
    publishDir "$params.OUTDIR/parameters", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR/tmle_inputs", mode: 'symlink', pattern: "*.csv"
    label "bigmem"
    
    input:
        path dataset
        path results

    output:
        path "permutation_dataset.arrow", emit: dataset
        path "*.yaml", emit: parameters

    script:
        limit = params.MAX_PERMUTATION_TESTS == null ? "" : "--limit=${params.MAX_PERMUTATION_TESTS}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /NegativeControl/bin/generate_permutation_data.jl \
        ${dataset} ${results} \
        ${limit} \
        --pval-col=${params.PVAL_COL} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --orders=${params.PERMUTATION_ORDERS} \
        --chunksize=${params.PERMUTATION_CHUNKSIZE} \
        --rng=${params.PERMUTATION_RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

process TMLE {
    container "olivierlabayle/targeted-estimation:up_tmle_dep"
    publishDir "$params.OUTDIR/csvs",  mode: 'symlink', pattern: "*.csv"
    publishDir "$params.OUTDIR/hdf5files/inf_curves",  mode: 'symlink', pattern: "*.hdf5"
    label "bigmem"
    label "multithreaded"

    input:
        path data
        path parameterfile
        path estimatorfile
    
    output:
        path "${outprefix}.csv", emit: tmle_csv
        path "${outprefix}.hdf5", optional: true, emit: inf_curve
    
    script:
        save_ic = params.NB_VAR_ESTIMATORS !== 0 ? '--save-ic' : ''
        outprefix = "tmle." + parameterfile.getName().replace(".yaml", "")
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --threads=${task.cpus} --startup-file=no /TargetedEstimation.jl/scripts/tmle.jl \
        $data $parameterfile $estimatorfile $outprefix \
        $save_ic \
        --pval-threshold=${params.PVAL_SIEVE}
        """
}


workflow {
    results_file = Channel.value(file(params.RESULTS_FILE))

    dataset, param_files = GeneratePermutationTestsData(results_file)
    RunPermutationTests(dataset, param_files)


}