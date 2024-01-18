using FeatureSetEnrichments
using VKCComputing
using EEGMicrobiome
using DataFrames
using CSV
using HypothesisTests
using MultipleTesting
using GLM

#kload data
concurrent_3m = load_cohort("concurrent_3m")
concurrent_6m = load_cohort("concurrent_6m")
concurrent_12m = load_cohort("concurrent_12m")
future_3m6m = load_cohort("future_3m6m")
future_3m12m = load_cohort("future_3m12m")
future_6m12m = load_cohort("future_6m12m")

concurrent_3m_func = DataFrame(get(concurrent_3m))

load_taxonomic_profiles!(mbo)
load_functional_profiles!(mbo)
unique!(mbo, :seqprep)
sort!(mbo, "age")
subset!(mbo, "age"=> ByRow(!ismissing))

tps = ("3m", "6m", "12m")
mbotps = Tuple((select(
            unique(
                subset(mbo, "visit"=> ByRow(v-> !ismissing(v) && v == "$(tp)o"), "subject"=> ByRow(!ismissing)),
                "subject"
            ), 
            "subject", "age"=>"stool_age", "seqprep"=>"sample", "age", Cols(:))
        for tp in tps
))

eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency"=>"amp")]

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)


