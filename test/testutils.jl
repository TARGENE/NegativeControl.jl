using TMLE
using TargetedEstimation

function make_estimates()
    IATE₁ = TMLE.TMLEstimate(
        estimand = IATE(
            outcome = "High light scatter reticulocyte percentage",
            treatment_values = (
                rs10043934 = (case="GA", control="GG"), 
                RSID_103 = (case="GA", control="GG")
            ),
            treatment_confounders = (
                rs10043934 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6), 
                RSID_103 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
            ),
            outcome_extra_covariates = ("Age-Assessment", "Genetic-Sex")
        ),
        estimate = -1.,
        std = 0.003,
        n = 10,
        IC = []
    )
    IATE₂ = TMLE.TMLEstimate(
        estimand = IATE(
            outcome = "High light scatter reticulocyte percentage",
            treatment_values = (
                rs10043934 = (case="GA", control="GG"), 
                RSID_103 = (case="AA", control="GA")
            ),
            treatment_confounders = (
                rs10043934 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6), 
                RSID_103 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
            ),
            outcome_extra_covariates = ("Age-Assessment", "Genetic-Sex")
        ),
        estimate = -0.003,
        std = 0.003,
        n = 10,
        IC = []
    )

    jointIATE = TMLE.ComposedEstimate(
        estimand = TMLE.ComposedEstimand(TMLE.joint_estimand, (IATE₁.estimand, IATE₂.estimand)),
        estimates = (IATE₁, IATE₂),
        estimate = [-1., -0.003],
        cov = [
        0.003 0.
        0. 0.003
        ],
        n = 10
    )

    ATE₁ = TMLE.TMLEstimate(
        estimand = ATE(
            outcome = "L50-L54 Urticaria and erythema",
            treatment_values = (
                rs117913124 = (case="GA", control="GG"), 
                RSID_104 = (case="GA", control="GG")
            ),
            treatment_confounders = (
                rs117913124 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6), 
                RSID_104 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
            ),
            outcome_extra_covariates = ("Age-Assessment", "Genetic-Sex")
        ),
        estimate = 0.003,
        std = 0.003,
        n = 20,
        IC = []
    )

    return [(TMLE=IATE₁,), (TMLE=IATE₂,), (TMLE=jointIATE,), (TMLE=ATE₁,)]
end

function save(estimates; prefix="tmle_output")
    outputs = TargetedEstimation.Outputs(
        json=TargetedEstimation.JSONOutput(filename=prefix*".json"),
        jls=TargetedEstimation.JLSOutput(filename=prefix*".jls"),
        hdf5=TargetedEstimation.HDF5Output(filename=prefix*".hdf5")
    )
    TargetedEstimation.initialize(outputs)
    batches = collect(Iterators.partition(estimates, 2))
    nbatches = length(batches)
    for (batchid, batch) in enumerate(batches)
        # Append JSON Output
        TargetedEstimation.update_file(outputs.json, batch; finalize=nbatches==batchid)
        # Append JLS Output
        TargetedEstimation.update_file(outputs.jls, batch)
        # Append HDF5 Output
        TargetedEstimation.update_file(outputs.hdf5, batch)
    end
end

make_fake_outputs(estimates_generator=make_estimates; prefix="tmle_output") = save(estimates_generator(); prefix=prefix)

function clean(;prefix="tmle_output")
    dir_, prefix_ = splitdir(prefix)
    dir = dir_ == "" ? "." : dir_
    for filename in readdir(dir)
        if startswith(filename, prefix_)
            rm(joinpath(dir_, filename))
        end
    end
end