using VKCComputing
using EEGMicrobiome
using FeatureSetEnrichments
using DataFrames
using CSV
using CairoMakie
using Preferences
using ColorSchemes

# y axis = gene abundances
# x axis = eeg feature
# scatter of samples, colored by bugs

fsea_results = CSV.read("data/outputs/fsea/true_ages_fsea.csv", DataFrame)

eeg = load_eeg()
eeg.peak_latency_P1_corrected = eeg.peak_latency_P1 .- eeg.peak_latency_N1
eeg.peak_latency_N2_corrected = eeg.peak_latency_N2 .- eeg.peak_latency_P1
eeg.peak_amp_P1_corrected = eeg.peak_amp_P1 .- eeg.peak_amp_N1
eeg.peak_amp_N2_corrected = eeg.peak_amp_N2 .- eeg.peak_amp_P1
rename!(eeg, "age" => "eeg_age")
mbo = load_microbiome(eeg.subject)
transform!(mbo, "visit" => ByRow(v -> ismissing(v) ? v : replace(v, "mo" => "m")) => "visit")

unique!(mbo, :seqprep)
sort!(mbo, "age")
subset!(mbo, "age" => ByRow(!ismissing))

tps = ("3m", "6m", "12m")
mbotps = Tuple((select(
    unique(
        subset(mbo, "visit" => ByRow(v -> !ismissing(v) && v == "$(tp)"), "subject" => ByRow(!ismissing)),
        "subject"
    ),
    "subject", "age" => "stool_age", "seqprep" => "sample", "age", Cols(:))
                for tp in tps
))

eeg_features = names(eeg, r"peak_")

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)

na_unirefs = Set(reduce(union, values(na_map)))
na_unirefs_full = Set(reduce(union, values(na_map_full)))

humann_files = mapreduce(vcat, readdir(joinpath(load_preference(VKCComputing, "mgx_analysis_dir"), "humann", "main"); join=true)) do f
    m = match(r"(SEQ\d+)", f)
    isnothing(m) && return DataFrame()
    sample = replace(basename(f), r"(SEQ\d+)_S\d+.+" => s"\1")
    sample ∈ mbo.seqprep || return DataFrame()
    if contains(basename(f), "genefamilies.tsv") && m[1] ∈ mbo.seqprep
        @info basename(f)
        df = CSV.read(f, DataFrame)
        rename!(df, ["feature", "abundance"])
        subset!(df, "feature" => ByRow(f -> contains(f, "|")))
        transform!(df, "feature" => ByRow(f -> begin
            (uniref, species) = split(f, "|")
            uniref = replace(uniref, "UniRef90_" => "")
            return (; uniref, species)
        end) => ["uniref", "species"])
        subset!(df, "uniref" => ByRow(u -> u ∈ na_unirefs))
        df.sample .= sample
        return df
    else
        return DataFrame()
    end
end

transform!(humann_files, "species" => ByRow(s -> begin
    s == "unclassified" && return (; genus="unclassified", species="unclassified")
    (g, s) = split(s, ".")

    contains(g, "_unclassified") && return (; genus=replace(g, "_unclassified"=>""), species=s)
    (; genus=g, species=s)
end) => ["genus", "species"]
)


#-

topgenera = mapreduce(vcat, eachrow(subset(fsea_results, "q₀" => ByRow(<(0.2))))) do row
    tp = row.timepoint
    gs = row.geneset
    eeg_feat = row.eeg_feature

    unirefs = na_map[gs]

    df = subset(humann_files, "uniref" => ByRow(u -> u ∈ unirefs))
    df = grouptop(df, 10)
    df = sort(combine(groupby(df, "genus"), "abundance"=>sum => "abundance"), "abundance"; rev=true)

    df.timepoint .= tp
    df.geneset .= gs
    df.eeg_feature .= eeg_feat
    df
end

topspecies = mapreduce(vcat, eachrow(subset(fsea_results, "q₀" => ByRow(<(0.2))))) do row
    tp = row.timepoint
    gs = row.geneset
    eeg_feat = row.eeg_feature

    unirefs = na_map[gs]

    df = subset(humann_files, "uniref" => ByRow(u -> u ∈ unirefs))
    df = grouptop(df, 10; groupcol="species")
    df = sort(combine(groupby(df, "species"), "abundance"=>sum => "abundance"), "abundance"; rev=true)

    df.timepoint .= tp
    df.geneset .= gs
    df.eeg_feature .= eeg_feat
    df
end

CSV.write("data/outputs/fsea_top_genera.csv", topgenera)
open("data/outputs/topspecies.txt", "w") do io
    println.(io, filter(!=("other"), unique(topspecies.species)))
end


CSV.write("data/outputs/fsea_top_species.csv", topspecies)
open("data/outputs/topgenera.txt", "w") do io
    println.(io, filter(!=("other"), unique(topgenera.genus)))
end

#-

for row in eachrow(subset(fsea_results, "q₀" => ByRow(<(0.2))))
    tp = row.timepoint
    gs = row.geneset
    eeg_feat = row.eeg_feature

    unirefs = na_map[gs]

    subdf = let df = subset(humann_files, "uniref" => ByRow(u -> u ∈ unirefs))
        df2 = leftjoin(mbotps[i], grouptop(df); on="sample")
        df2 = dropmissing(leftjoin!(select(df2, "visit"=> "timepoint", Cols(:)), select(eeg, "subject", "timepoint", eeg_feat); on=["subject", "timepoint"]))
    end

    cs = Dict(bug=> i for (i, bug) in enumerate(sort(unique(subdf."genus"))))
    subdf.color = map(x -> cs[x], subdf.genus)

    fig = Figure(; size=(800, 600));

    es = round(row.es, digits=4)
    q₀ = round(row.q₀, digits=4)

    ax = Axis(fig[1,1]; xlabel=eeg_feat, ylabel="$gs abundance (RPKM)", title="E.S. $es, q₀ = $q₀")
    scatter!(ax, subdf[!, eeg_feat], subdf[!, "abundance"]; color=subdf.color, colormap=(:tab10, 0.6))
    
    Legend(fig[1,2], [MarkerElement(; color=(ColorSchemes.tab10[i], 0.6), marker=:circle) for i in 1:length(keys(cs))],
                    sort(unique(subdf."genus"))
    )
    save(joinpath("data", "figures", "bugs", "$(tp)_$(gs)_$(eeg_feat).png"), fig)
end
    
#-

specs = Set(String[])
genera = Set(String[])
foreach(readdir(joinpath(load_preference(VKCComputing, "mgx_analysis_dir"), "humann", "main"); join=true)) do f
    m = match(r"(SEQ\d+)", f)
    isnothing(m) && return nothing
    sample = replace(basename(f), r"(SEQ\d+)_S\d+.+" => s"\1")
    sample ∈ mbo.seqprep || return nothing
    if contains(basename(f), "genefamilies.tsv") && m[1] ∈ mbo.seqprep
        @info basename(f)
        spec = CSV.read(f, DataFrame)[!, 1]
        filter!(spec -> contains(spec, "|"), spec)
        newspecs = Set(split(s, "|")[2] for s in spec)
        union!(specs, newspecs)
        union!(genera, Set(split(s, ".")[1] for s in newspecs))
    end
end
