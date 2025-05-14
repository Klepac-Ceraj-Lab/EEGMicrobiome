using VKCComputing
using EEGMicrobiome
using FeatureSetEnrichments
using DataFrames
using CSV
using Chain
using CairoMakie
using Preferences
using ColorSchemes

# y axis = gene abundances
# x axis = eeg feature
# scatter of samples, colored by bugs

fsea_results = CSV.read("data/outputs/fsea/true_ages_fsea.csv", DataFrame)

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)

na_unirefs = Set(reduce(union, values(na_map)))
na_unirefs_full = Set(reduce(union, values(na_map_full)))

tps = ("v1", "v2", "v3")
ftps = ("v1v2", "v1v3", "v2v3")

mdata = load_cohorts()
mdata.taxprofile = map(mdata.taxprofile) do f
    ismissing(f) && return missing
    joinpath(dirname(f), "mpa_v31_CHOCOPhlAn_201901",
        replace(basename(f), r"_profile\.tsv" => "_mpa_v31_CHOCOPhlAn_201901_profile.tsv")
        )
end

long_sub = let
  wide_sub = select(
    leftjoin(
      select(unstack(mdata, "subject_id", "visit", "eeg_age"),
        "subject_id", "v1" => "eeg_v1", "v2" => "eeg_v2", "v3" => "eeg_v3"),
      select(unstack(mdata, "subject_id", "visit", "stool_age"),
        "subject_id", "v1" => "seqprep_v1", "v2" => "seqprep_v2", "v3" => "seqprep_v3"),
      on="subject_id"),
    "subject_id", r"v1", r"v2", r"v3"
  )

  long_sub = DataFrame()
  for row in eachrow(wide_sub), tp in tps
    stool_age = row["seqprep_$tp"]
    eeg_age = row["eeg_$tp"]
    push!(long_sub, (; subject_id=row.subject_id, timepoint=tp, stool_age, eeg_age); cols=:union)
  end

  @chain long_sub begin
    subset!(AsTable(["stool_age", "eeg_age"]) => ByRow(nt -> !all(ismissing, nt)))
    transform!(AsTable(["stool_age", "eeg_age"]) => ByRow(nt -> minimum(skipmissing(values(nt)))) => "minage")
    sort!("minage")
  end
end


v1 = get_cohort(mdata, "v1")
v2 = get_cohort(mdata, "v2")
v3 = get_cohort(mdata, "v3")
v1v2 = get_cohort(mdata, "v1v2")
v1v3 = get_cohort(mdata, "v1v3")
v2v3 = get_cohort(mdata, "v2v3")

humann_files = let allseqs = Set(mapreduce(df-> df.seqprep, vcat, (v1,v2,v3, v1v2, v1v3, v2v3)))
    mapreduce(vcat, readdir(joinpath(load_preference(VKCComputing, "mgx_analysis_dir"), "humann", "main"); join=true)) do f
        m = match(r"(SEQ\d+)", f)
        isnothing(m) && return DataFrame()
        sample = replace(basename(f), r"(SEQ\d+)_S\d+.+" => s"\1")
        sample ∈ allseqs || return DataFrame()
        if contains(basename(f), "genefamilies.tsv") && m[1] ∈ allseqs
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
end;


CSV.write("data/outputs/fsea_top_species.csv", topspecies)
open("data/outputs/topgenera.txt", "w") do io
    println.(io, filter(!=("other"), unique(topgenera.genus)))
end;

colormap_age = :viridis
colors_timepoints = Dict(tp => c for (tp, c) in zip(("3m", "6m", "12m"), cgrad(colormap_age)[[0.0, 0.4, 0.8]]))
#-
for grp in groupby(subset(fsea_results, "q₀" => ByRow(<(0.2))), ["eeg_feature"])
    fig = Figure(; size=(1600,1200))
    sqiter = ceil(Int, sqrt(nrow(grp)))

    for (row, (i,j)) in zip(eachrow(sort(grp, "timepoint"; lt= (x,y) -> x == "3m" || y == "12m")), Iterators.product(1:sqiter, 1:sqiter))
        grid = GridLayout(fig[j,i])
        tp = row.timepoint
        gs = row.geneset
        eeg_feat = row.eeg_feature

        unirefs = na_map[gs]

        subdf = let df = subset(humann_files, "uniref" => ByRow(u -> u ∈ unirefs))
            df2 = leftjoin(vcat(v1,v2,v3), grouptop(df); on="seqprep"=>"sample")
            dropmissing!(df2)
        end

        cs = Dict(bug=> i for (i, bug) in enumerate(sort(unique(subdf."genus"))))
        subdf.color = map(x -> cs[x], subdf.genus)

        es = round(row.es, digits=4)
        q₀ = round(row.q₀, digits=4)

        local ax = Axis(grid[1,1]; xlabel=eeg_feat, ylabel="$gs abundance (RPKM)", title="E.S. $es, q₀ = $q₀", titlecolor=colors_timepoints[row.timepoint])
        scatter!(ax, subdf[!, eeg_feat], subdf[!, "abundance"]; color=subdf.color, colormap=(:tab10, 0.6))
        
        Legend(grid[1, 2], [MarkerElement(; color=(ColorSchemes.tab10[i], 0.6), marker=:circle) for i in 1:length(keys(cs))],
                        sort(unique(subdf."genus"))
        )
    end

    save(joinpath("data", "figures", "bugs", "$(first(grp.eeg_feature)).png"), fig)
    save(joinpath("data", "figures", "bugs", "$(first(grp.eeg_feature)).svg"), fig)
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
