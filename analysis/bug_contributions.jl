using VKCComputing
using EEGMicrobiome
using FeatureSetEnrichments
using DataFrames
using CSV
using CairoMakie
using Preferences

# y axis = gene abundances
# x axis = eeg feature
# scatter of samples, colored by bugs

eeg = load_eeg()
eeg.peak_latency_P1_corrected = eeg.peak_latency_P1 .- eeg.peak_latency_N1
eeg.peak_latency_N2_corrected = eeg.peak_latency_N2 .- eeg.peak_latency_P1
eeg.peak_amp_P1_corrected = eeg.peak_amp_P1 .- eeg.peak_amp_N1
eeg.peak_amp_N2_corrected = eeg.peak_amp_N2 .- eeg.peak_amp_P1
rename!(eeg, "age"=> "eeg_age")
mbo = load_microbiome(eeg.subject)

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

eeg_features = names(eeg, r"peak_")

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)

humann_files = Dict(seq => df for (seq, df) in filter(!isnothing, 
    map(readdir(joinpath(load_preference(VKCComputing, "mgx_analysis_dir"), "humann", "main"); join=true)) do f
        m = match(r"(SEQ\d+)", f)
        isnothing(m) && return nothing

        if contains(basename(f), "genefamilies.tsv") && m[1] ∈ mbo.seqprep
            return (String(m[1]), CSV.read(f, DataFrame))
        else
            return nothing
        end
end))

#-

eeg_feat = "peak_amp_N1"
gs = "Glutamate synthesis"
tp = "3m"

mbotp = mbotps[1]
mbotp.seqprep

