using FeatureSetEnrichments
using VKCComputing
using EEGMicrobiome
using DataFrames
using CSV

# load data
eeg = load_eeg()
eeg.peak_latency_P1_corrected = eeg.peak_latency_P1 .- eeg.peak_latency_N1
eeg.peak_latency_N2_corrected = eeg.peak_latency_N2 .- eeg.peak_latency_P1
eeg.peak_amp_P1_corrected = eeg.peak_amp_P1 .- eeg.peak_amp_N1
eeg.peak_amp_N2_corrected = eeg.peak_amp_N2 .- eegvep.peak_amp_P1

mbo = load_microbiome(eeg.subject)

load_taxonomic_profiles!(mbo)
load_functional_profiles!(mbo)
unique!(mbo, :seqprep)
sort!(mbo, "age")
subset!(mbo, "age"=> ByRow(!ismissing))

tps = ("3m", "6m", "12m")
eeg_features = names(eeg, r"peak_")

#-

for tp in tps
    tpi = parse(Int, match(r"\d+", tp).match)
    subeeg = subset(eeg, "timepoint"=> ByRow(==(tp)))
    leftjoin!(subeeg, select(unique(subset(mbo, "visit"=> ByRow(v-> !ismissing(v) && v == "$(tp)o"), "subject"=> ByRow(!ismissing)), :subject), "subject", "age"=>"stool_age", "seqprep"=>"sample", Cols(Not("age"))); on=:subject)
    subset!(subeeg, "seqprep"=> ByRow(!ismissing))
    for feature in eeg_features
        EEGMicrobiome.runlms(subeeg, "./data/outputs/lms/$(feature)_$(tp)_lms.csv", feature, names(subeeg, r"^UniRef"))
    end
end

## Run FSEA

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)

fsea_df = let fsea_df = DataFrame()
    for tp in tps, feature in eeg_features
        filepath = "./data/outputs/lms/$(feature)_$(tp)_lms.csv"
        lms = CSV.read(filepath, DataFrame)

        for (key, unirefs) in pairs(na_map)
            s = Set("UniRef90_$uniref" for uniref in unirefs)
            lms[!, key] = lms.feature .∈ Ref(s)
        end

        for (key, unirefs) in pairs(na_map)
            idx = findall(lms[!, key])
            sum(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.Permutation(5000), lms.z, idx)
            push!(fsea_df, (; timepoint=tp, eeg_feature=feature, geneset=key, pvalue=pvalue(result), es = enrichment_score(result)))
        end
    end
    fsea_df
end

transform!(fsea_df, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(fsea_df, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(fsea_df, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")

CSV.write("data/outputs/fsea/true_ages_fsea.csv", fsea_df)