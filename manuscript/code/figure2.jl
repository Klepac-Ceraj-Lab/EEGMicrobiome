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
load_functional_profiles!(concurrent_3m_func);

concurrent_6m_func = DataFrame(get(concurrent_6m))
load_functional_profiles!(concurrent_6m_func);

concurrent_12m_func = DataFrame(get(concurrent_12m))
load_functional_profiles!(concurrent_12m_func);

future_3m6m_func = DataFrame(get(future_3m6m))
transform!(future_3m6m_func, AsTable(["age", "eeg_age"])=> ByRow(nt-> nt.eeg_age - nt.age)=> "age_diff")
load_functional_profiles!(future_3m6m_func);

future_3m12m_func = DataFrame(get(future_3m12m))
load_functional_profiles!(future_3m12m_func);

future_6m12m_func = DataFrame(get(future_6m12m))
load_functional_profiles!(future_6m12m_func);

##


eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency"=>"amp")]

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)


for feature in eeg_features
    EEGMicrobiome.runlms(concurrent_3m_func, "./data/outputs/lms/$(feature)_3m_lms.csv",
                         feature, names(concurrent_3m_func, r"^UniRef")
    )
end

for feature in eeg_features
    EEGMicrobiome.runlms(concurrent_6m_func, "./data/outputs/lms/$(feature)_6m_lms.csv",
                         feature, names(concurrent_6m_func, r"^UniRef")
    )
end

for feature in eeg_features
    EEGMicrobiome.runlms(concurrent_12m_func, "./data/outputs/lms/$(feature)_12m_lms.csv",
                         feature, names(concurrent_12m_func, r"^UniRef")
    )
end

