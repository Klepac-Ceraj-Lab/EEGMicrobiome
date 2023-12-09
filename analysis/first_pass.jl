using FeatureSetEnrichments
using VKCComputing
using EEGMicrobiome
using DataFrames
using CSV
using HypothesisTests
using MultipleTesting
using GLM

# load data
eeg = load_eeg()
eeg.peak_latency_P1_corrected = eeg.peak_latency_P1 .- eeg.peak_latency_N1
eeg.peak_latency_N2_corrected = eeg.peak_latency_N2 .- eeg.peak_latency_P1
eeg.peak_amp_P1_corrected = eeg.peak_amp_P1 .- eeg.peak_amp_N1
eeg.peak_amp_N2_corrected = eeg.peak_amp_N2 .- eeg.peak_amp_P1
rename!(eeg, "age"=> "eeg_age")
mbo = load_microbiome(eeg.subject)

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

eeg_features = names(eeg, r"peak_")

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)

## Run Linear models

### For concurrent timepoints


for (tp, df) in zip(tps, mbotps)
    tpi = parse(Int, match(r"\d+", tp).match)
    subeeg = subset(eeg, "timepoint"=> ByRow(==(tp)))
    leftjoin!(subeeg, df; on=:subject)
    subset!(subeeg, "sample"=> ByRow(!ismissing))
    for feature in eeg_features
        EEGMicrobiome.runlms(subeeg, "./data/outputs/lms/$(feature)_$(tp)_lms.csv", feature, names(subeeg, r"^UniRef"))
    end
end

### for Future prediction

#### 3m Microbiome -> 6m EEG

let
    tp = "3m"
    df = mbotps[1]
    subeeg = subset(eeg, "timepoint"=> ByRow(==("6m")))
    leftjoin!(subeeg, df; on=:subject)
    transform!(subeeg, AsTable(["age", "eeg_age"])=> ByRow(nt-> nt.eeg_age - nt.age) => "age_diff")
    subset!(subeeg, "sample"=> ByRow(!ismissing))
    for feature in eeg_features
        EEGMicrobiome.runlms(subeeg, "./data/outputs/lms/$(feature)_$(tp)_future6m_lms.csv", feature, names(subeeg, r"^UniRef");
            additional_cols = [:age_diff],
            formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
        )
    end
end

#### 3m Microbiome -> 12m EEG

let
    tp = "3m"
    df = mbotps[1]
    subeeg = subset(eeg, "timepoint"=> ByRow(==("12m")))
    leftjoin!(subeeg, df; on=:subject)
    subset!(subeeg, "sample"=> ByRow(!ismissing))
    transform!(subeeg, AsTable(["age", "eeg_age"])=> ByRow(nt-> nt.eeg_age - nt.age) => "age_diff")
    for feature in eeg_features
        EEGMicrobiome.runlms(subeeg, "./data/outputs/lms/$(feature)_$(tp)_future12m_lms.csv", feature, names(subeeg, r"^UniRef");
            additional_cols = [:age_diff],
            formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
        )
    end
end


#### 6m Microbiome -> 12m EEG

let
    tp = "6m"
    df = mbotps[2]
    subeeg = subset(eeg, "timepoint"=> ByRow(==("12m")))
    leftjoin!(subeeg, df; on=:subject)
    subset!(subeeg, "sample"=> ByRow(!ismissing))
    transform!(subeeg, AsTable(["age", "eeg_age"])=> ByRow(nt-> nt.eeg_age - nt.age) => "age_diff")
    for feature in eeg_features
        EEGMicrobiome.runlms(subeeg, "./data/outputs/lms/$(feature)_$(tp)_future12m_lms.csv", feature, names(subeeg, r"^UniRef");
            additional_cols = [:age_diff],
            formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
        )
    end
end

## Run FSEA


#-

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

#- 

## Run FSEA

future12m_fsea_df = let fsea_df = DataFrame()
    for tp in ("3m", "6m"), feature in eeg_features
        # contains(feature, "corrected") && continue
        filepath = "./data/outputs/lms/$(feature)_$(tp)_future12m_lms.csv"
        lms = CSV.read(filepath, DataFrame)

        for (key, unirefs) in pairs(na_map)
            s = Set("UniRef90_$uniref" for uniref in unirefs)
            lms[!, key] = lms.feature .∈ Ref(s)
        end

        for (key, unirefs) in pairs(na_map)
            idx = findall(lms[!, key])
            sum(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.Permutation(5000), lms.z, idx)
            push!(fsea_df, (; timepoint="$(tp)_future12m", eeg_feature=feature, geneset=key, pvalue=pvalue(result), es = enrichment_score(result)))
        end
    end
    fsea_df
end

transform!(future12m_fsea_df, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(future12m_fsea_df, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(future12m_fsea_df, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")

CSV.write("data/outputs/fsea/future12m_true_ages_fsea.csv", future12m_fsea_df)
# future12m_fsea_df = CSV.read("data/outputs/fsea/future12m_true_ages_fsea.csv", DataFrame)

#-

future6m_fsea_df = let fsea_df = DataFrame()
    for feature in eeg_features
        # contains(feature, "corrected") && continue
        tp = "3m"
        filepath = "./data/outputs/lms/$(feature)_$(tp)_future6m_lms.csv"
        lms = CSV.read(filepath, DataFrame)

        for (key, unirefs) in pairs(na_map)
            s = Set("UniRef90_$uniref" for uniref in unirefs)
            lms[!, key] = lms.feature .∈ Ref(s)
        end

        for (key, unirefs) in pairs(na_map)
            idx = findall(lms[!, key])
            sum(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.Permutation(5000), lms.z, idx)
            push!(fsea_df, (; timepoint="$(tp)_future6m", eeg_feature=feature, geneset=key, pvalue=pvalue(result), es = enrichment_score(result)))
        end
    end
    fsea_df
end

transform!(future6m_fsea_df, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(future6m_fsea_df, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(future6m_fsea_df, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")

CSV.write("data/outputs/fsea/future6m_true_ages_fsea.csv", future6m_fsea_df)
# future6m_fsea_df = CSV.read("data/outputs/fsea/future6m_true_ages_fsea.csv", DataFrame)

#- 

## Visualize FSEA

using CairoMakie

unirefs_3m = replace.(Iterators.filter(
                        n -> sum(skipmissing(mbotps[1][!,n])) > 0, 
                        names(mbotps[1], r"^UniRef")),
                "UniRef90_" => ""
)
unirefs_6m = replace.(Iterators.filter(
                        n -> sum(skipmissing(mbotps[2][!,n])) > 0, 
                        names(mbotps[1], r"^UniRef")),
                "UniRef90_" => ""
)
unirefs_12m = replace.(Iterators.filter(
                        n -> sum(skipmissing(mbotps[3][!,n])) > 0, 
                        names(mbotps[1], r"^UniRef")),
                "UniRef90_" => ""
)
values_full = values(na_map_full)
values_consolidate = values(na_map)

#-

fig = Figure(; size = (900, 900))
ax1 = Axis(fig[1,1]; title = "genes in gene sets", xlabel = "length", ylabel = "frequency")
ax2 = Axis(fig[2,1]; title = "genes in (consolidated) gene sets", xlabel = "length", ylabel = "frequency")
# histogram of number of genes per geneset


density!(ax1, Vector{Int}(filter(>(0), map(value -> length(intersect(value, unirefs_3m)), na_map_full)) |> collect))
density!(ax2, Vector{Int}(filter(>(0), map(value -> length(intersect(value, unirefs_3m)), na_map)) |> collect))
density!(ax1, Vector{Int}(filter(>(0), map(value -> length(intersect(value, unirefs_6m)), na_map_full)) |> collect))
density!(ax2, Vector{Int}(filter(>(0), map(value -> length(intersect(value, unirefs_6m)), na_map)) |> collect))
density!(ax1, Vector{Int}(filter(>(0), map(value -> length(intersect(value, unirefs_12m)), na_map_full)) |> collect))
density!(ax2, Vector{Int}(filter(>(0), map(value -> length(intersect(value, unirefs_12m)), na_map)) |> collect))

Legend(fig[3,1:2], [MarkerElement(; color=c, marker=:rect) for c in Makie.wong_colors()[1:3]],
                 ["3m", "6m", "12m"], "Genes from...";
                 tellwidth=false, tellheight=true, orientation=:horizontal
)

ax3 = Axis(fig[1,2]; xticks=(1:3, [tps...]), ylabel="N samples")
barplot!(ax3, 1:3, [nrow(df) for df in mbotps]; color=Makie.wong_colors()[1:3])

ax4 = Axis(fig[2,2]; xlabel="peak amplitude N1", ylabel="peak latency N1")

for (i, tp) in enumerate(tps)
    df = subset(select(eeg, first(eeg_features, 2)..., "timepoint"), "timepoint"=> ByRow(==(tp)))
    scatter!(ax4, df[!, 1], df[!, 2]; color = (Makie.wong_colors()[i], 0.6))
end

save("./data/figures/fsea_diagnostics.png", fig)
fig


#- 

fig = Figure(; size = (900, 900))

for i in eachindex(eeg_features), j in eachindex(eeg_features)
    i < j || continue
    ax = Axis(fig[i,j]; xlabel=eeg_features[i], ylabel=eeg_features[j])
    for (i, tp) in enumerate(tps)
        df = subset(select(eeg, eeg_features[i], eeg_features[j], "timepoint"), "timepoint"=> ByRow(==(tp)))
        scatter!(ax4, df[!, 1], df[!, 2]; color = (Makie.wong_colors()[i], 0.6))
    end
end

fig

#-

let
    tp = "3m"
    df = mbotps[1]
    subeeg = subset(eeg, "timepoint"=> ByRow(==("6m")))
    leftjoin!(subeeg, df; on=:subject)
    subset!(subeeg, "sample"=> ByRow(!ismissing))
    @info "$tp microbiome vs 6m eeg, N = $(size(subeeg, 1))"
end
let
    tp = "3m"
    df = mbotps[1]
    subeeg = subset(eeg, "timepoint"=> ByRow(==("12m")))
    leftjoin!(subeeg, df; on=:subject)
    subset!(subeeg, "sample"=> ByRow(!ismissing))
    @info "$tp microbiome vs 12m eeg, N = $(size(subeeg, 1))"
end
let
    tp = "6m"
    df = mbotps[2]
    subeeg = subset(eeg, "timepoint"=> ByRow(==("12m")))
    leftjoin!(subeeg, df; on=:subject)
    subset!(subeeg, "sample"=> ByRow(!ismissing))
    @info "$tp microbiome vs 12m eeg, N = $(size(subeeg, 1))"
end