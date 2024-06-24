using FeatureSetEnrichments
using VKCComputing
using EEGMicrobiome
using Microbiome
using DataFrames
using Distances
using CSV
using HypothesisTests
using MultipleTesting
using Clustering
using GLM
using Distributions
using CairoMakie
using ColorSchemes
using ThreadsX
using SparseArrays

mdata = load_cohorts()

v1 = get_cohort(mdata, "v1")
v2 = get_cohort(mdata, "v2")
v3 = get_cohort(mdata, "v3")
v1v2 = get_cohort(mdata, "v1v2")
v1v3 = get_cohort(mdata, "v1v3")
v2v3 = get_cohort(mdata, "v2v3")

unirefs_by_sample = let 
	files = String.(skipmissing(mdata.genefamilies))
	mapreduce(vcat, files) do f
		df = CSV.read(f, DataFrame)
		rename!(df, ["feature", "abundance"])
		subset!(df, "feature"=> ByRow(f-> !contains(f, '|')))
		sample = replace(basename(f), r"(SEQ\d+).+"=> s"\1")
		df.sample .= sample
		df
    end
end

unirefs = let 
    fs = unique(unirefs_by_sample.feature)
    ss = unique(unirefs_by_sample.sample)
    fsmap = Dict(f=> i for (i, f) in enumerate(fs))
    ssmap = Dict(s=> i for (i, s) in enumerate(ss))

    mat = spzeros(length(fs), length(ss))
    foreach(eachrow(unirefs_by_sample)) do row
	mat[fsmap[row.feature], ssmap[row.sample]] = row.abundance
    end
    CommunityProfile(mat, GeneFunction.(fs), MicrobiomeSample.(ss))
end


set!(unirefs, select(mdata, "seqprep"=> "sample", Cols(:)))


v1_func = let comm = unirefs[:, v1.seqprep]
    mat = collect(abundances(comm)')
    ixs = ThreadsX.map(col-> sum(col) > 0, eachcol(mat))
    hcat(v1,  DataFrame(mat[:, ixs], featurenames(comm)[ixs]))
end

v2_func = let comm = unirefs[:, v2.seqprep]
    mat = collect(abundances(comm)')
    ixs = ThreadsX.map(col-> sum(col) > 0, eachcol(mat))
    hcat(v2,  DataFrame(mat[:, ixs], featurenames(comm)[ixs]))
end


v3_func = let comm = unirefs[:, v3.seqprep]
    mat = collect(abundances(comm)')
    ixs = ThreadsX.map(col-> sum(col) > 0, eachcol(mat))
    hcat(v3,  DataFrame(mat[:, ixs], featurenames(comm)[ixs]))
end

v1v2_func = let comm = unirefs[:, v1v2.seqprep]
    mat = collect(abundances(comm)')
    ixs = ThreadsX.map(col-> sum(col) > 0, eachcol(mat))
    hcat(v1v2, DataFrame(mat[:, ixs], featurenames(comm)[ixs]))
end

v1v3_func = let comm = unirefs[:, v1v3.seqprep]
    mat = collect(abundances(comm)')
    ixs = ThreadsX.map(col-> sum(col) > 0, eachcol(mat))
    hcat(v1v3, DataFrame(mat[:, ixs], featurenames(comm)[ixs]))
end

v2v3_func = let comm = unirefs[:, v2v3.seqprep]
    mat = collect(abundances(comm)')
    ixs = ThreadsX.map(col-> sum(col) > 0, eachcol(mat))
    hcat(v2v3, DataFrame(mat[:, ixs], featurenames(comm)[ixs]))
end


##


eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency"=>"amp")]

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v1_func, "./data/outputs/lms/$(feature)_nodiff_v1_lms.csv",
                         feature, names(v1_func, r"^UniRef"); age_col="stool_age"
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v2_func, "./data/outputs/lms/$(feature)_nodiff_v2_lms.csv",
                         feature, names(v2_func, r"^UniRef"); age_col="stool_age"
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v3_func, "./data/outputs/lms/$(feature)_nodiff_v3_lms.csv",
                         feature, names(v3_func, r"^UniRef"); age_col="stool_age"
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v1_func, "./data/outputs/lms/$(feature)_v1_lms.csv",
        feature, names(v1_func, r"^UniRef"),
        additional_cols = [:age_diff], age_col = "stool_age",
        formula = term(:func) ~ term(feature) + term(:stool_age) + term(:age_diff) + term(:n_trials)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v2_func, "./data/outputs/lms/$(feature)_v2_lms.csv",
        feature, names(v2_func, r"^UniRef"),
        additional_cols = [:age_diff], age_col = "stool_age",
        formula = term(:func) ~ term(feature) + term(:stool_age) + term(:age_diff) + term(:n_trials)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v3_func, "./data/outputs/lms/$(feature)_v3_lms.csv",
        feature, names(v3_func, r"^UniRef"),
        additional_cols = [:age_diff], age_col = "stool_age",
        formula = term(:func) ~ term(feature) + term(:stool_age) + term(:age_diff) + term(:n_trials)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v1v2_func, "./data/outputs/lms/$(feature)_v1v2_lms.csv", feature, names(v1v2_func, r"^UniRef");
        additional_cols = [:age_diff], age_col = "stool_age",
        formula = term(:func) ~ term(feature) + term(:stool_age) + term(:age_diff) + term(:n_trials)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v1v3_func, "./data/outputs/lms/$(feature)_v1v3_lms.csv", feature, names(v1v3_func, r"^UniRef");
        additional_cols = [:age_diff], age_col = "stool_age",
        formula = term(:func) ~ term(feature) + term(:stool_age) + term(:age_diff) + term(:n_trials)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(v2v3_func, "./data/outputs/lms/$(feature)_v2v3_lms.csv", feature, names(v2v3_func, r"^UniRef");
        additional_cols = [:age_diff], age_col = "stool_age",
        formula = term(:func) ~ term(feature) + term(:stool_age) + term(:age_diff) + term(:n_trials)
    )
end

## Run FSEA

tps = ("v1", "v2", "v3")

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
            length(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.Permutation(5000), lms.z, idx)
            push!(fsea_df, (; timepoint=tp,
                              eeg_feature = feature,
                              geneset     = key,
                              pvalue      = pvalue(result),
                              es          = enrichment_score(result),
                              ranks       = idx,
                              nfeatures   = size(lms, 1)))
        end
    end
    fsea_df
end

transform!(fsea_df, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(fsea_df, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(fsea_df, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")
sort!(fsea_df, "q₀")
CSV.write("data/outputs/fsea/concurrent_consolidated_fsea.csv", fsea_df)
# fsea_df = CSV.read("data/outputs/fsea/concurrent_consolidated_fsea.csv", DataFrame)
let res = Dict()
    for tp in tps, feature in eeg_features
        @info tp
        filepath = "./data/outputs/lms/$(feature)_$(tp)_lms.csv"
        lms = CSV.read(filepath, DataFrame)

        for (key, unirefs) in pairs(na_map)
            s = Set("UniRef90_$uniref" for uniref in unirefs)
            lms[!, key] = lms.feature .∈ Ref(s)
        end

        for (key, unirefs) in pairs(na_map)
            idx = findall(lms[!, key])
            res[key] = get(res, key, Dict())
            res[key][tp] = length(idx)
        end
    end
    res
end

fsea_df_u = let fsea_df = DataFrame()
    for tp in tps, feature in eeg_features
        filepath = "./data/outputs/lms/$(feature)_$(tp)_lms.csv"
        lms = CSV.read(filepath, DataFrame)

        for (key, unirefs) in pairs(na_map)
            s = Set("UniRef90_$uniref" for uniref in unirefs)
            lms[!, key] = lms.feature .∈ Ref(s)
        end

        for (key, unirefs) in pairs(na_map)
            idx = findall(lms[!, key])
            length(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.MWU(), lms.z, idx)
            push!(fsea_df, (; timepoint=tp,
                              eeg_feature = feature,
                              geneset     = key,
                              pvalue      = pvalue(result),
                              es          = enrichment_score(result),
                              ranks       = idx,
                              nfeatures   = size(lms, 1)))
        end
    end
    fsea_df
end

transform!(fsea_df_u, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(fsea_df_u, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(fsea_df_u, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")
sort!(fsea_df_u, "q₀")
CSV.write("data/outputs/fsea/concurrent_consolidated_fsea_u.csv", fsea_df_u)
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
            length(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.Permutation(5000), lms.z, idx)
            push!(fsea_df, (; timepoint="$(tp)_future12m",
                              eeg_feature = feature,
                              geneset     = key,
                              pvalue      = pvalue(result),
                              es          = enrichment_score(result),
                              ranks       = idx,
                              nfeatures   = size(lms, 1)))
        end
    end
    fsea_df
end

transform!(future12m_fsea_df, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(future12m_fsea_df, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(future12m_fsea_df, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")
sort!(future12m_fsea_df, "q₀")
CSV.write("data/outputs/fsea/future12m_consolidated_fsea.csv", future12m_fsea_df)
# future12m_fsea_df = CSV.read("data/outputs/fsea/future12m_consolidated_fsea.csv", DataFrame)

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
            @warn key
            idx = findall(lms[!, key])
            length(idx) > 5 || continue
            @info "keeping"
            result = fsea(FeatureSetEnrichments.Permutation(5000), lms.z, idx)
            push!(fsea_df, (; timepoint="$(tp)_future6m",
                              eeg_feature = feature,
                              geneset     = key,
                              pvalue      = pvalue(result),
                              es          = enrichment_score(result),
                              ranks       = idx,
                              nfeatures   = size(lms, 1)))

         end
    end
    fsea_df
end

transform!(future6m_fsea_df, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(future6m_fsea_df, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(future6m_fsea_df, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")
sort!(future6m_fsea_df, "q₀")
CSV.write("data/outputs/fsea/future6m_consolidated_fsea_u.csv", future6m_fsea_df)
# future6m_fsea_df = CSV.read("data/outputs/fsea/future6m_consolidated_fsea_u.csv", DataFrame)
#
future12m_fsea_df_u = let fsea_df_u = DataFrame()
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
            length(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.MWU(), lms.z, idx)
            push!(fsea_df_u, (; timepoint="$(tp)_future12m",
                              eeg_feature = feature,
                              geneset     = key,
                              pvalue      = pvalue(result),
                              es          = enrichment_score(result),
                              ranks       = idx,
                              nfeatures   = size(lms, 1)))
        end
    end
    fsea_df_u
end

transform!(future12m_fsea_df_u, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(future12m_fsea_df_u, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(future12m_fsea_df_u, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")
sort!(future12m_fsea_df_u, "q₀")
CSV.write("data/outputs/fsea/future12m_consolidated_fsea_u.csv", future12m_fsea_df_u)
# future12m_fsea_df_u = CSV.read("data/outputs/fsea/future12m_consolidated_fsea_u.csv", DataFrame)

#-

future6m_fsea_df_u = let fsea_df_u = DataFrame()
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
            length(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.MWU(), lms.z, idx)
            push!(fsea_df_u, (; timepoint="$(tp)_future6m",
                              eeg_feature = feature,
                              geneset     = key,
                              pvalue      = pvalue(result),
                              es          = enrichment_score(result),
                              ranks       = idx,
                              nfeatures   = size(lms, 1)))

         end
    end
    fsea_df_u
end

transform!(future6m_fsea_df_u, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(future6m_fsea_df_u, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(future6m_fsea_df_u, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")
sort!(future6m_fsea_df_u, "q₀")
CSV.write("data/outputs/fsea/future6m_consolidated_fsea_u.csv", future6m_fsea_df_u)
# future6m_fsea_df_u = CSV.read("data/outputs/fsea/future6m_consolidated_fsea_u.csv", DataFrame)
############
# Plotting #
############

fsea_df = CSV.read("data/outputs/fsea/concurrent_consolidated_fsea.csv", DataFrame)
future6m_fsea_df = CSV.read("data/outputs/fsea/future6m_consolidated_fsea.csv", DataFrame)
future12m_fsea_df = CSV.read("data/outputs/fsea/future12m_consolidated_fsea.csv", DataFrame)

for row in eachrow(subset(fsea_df, "q₀" => ByRow(<(0.2))))
    res = fsea(1:row.nfeatures, row.ranks)
    fig = Figure(; size=(800, 800))
    gr = GridLayout(fig[1,1])
    _, ax1, ax2 = plot!(gr, res)

    text!(ax1, 0, 1;
            text = "$(row.timepoint)\n$(row.eeg_feature)\n$(row.geneset)",
            align = (:left, :top),
            offset = (4, -2),
            space = :relative,
            fontsize = 12
    )
    save("data/figures/enrichments/concurrent_$(row.timepoint)_$(row.eeg_feature)_$(replace(row.geneset, " "=>"-")).png", fig)
end

for row in eachrow(subset(future12m_fsea_df, "q₀" => ByRow(<(0.2))))
    res = fsea(1:row.nfeatures, row.ranks)
    fig = Figure(; size=(800, 800))
    gr = GridLayout(fig[1,1])
    _, ax1, ax2 = plot!(gr, res)

    text!(ax1, 0, 1;
            text = "$(row.timepoint)\n$(row.eeg_feature)\n$(row.geneset)",
            align = (:left, :top),
            offset = (4, -2),
            space = :relative,
            fontsize = 12
    )
    save("data/figures/enrichments/future_$(row.timepoint)_$(row.eeg_feature)_$(replace(row.geneset, " "=>"-")).png", fig)
end

for row in eachrow(subset(future6m_fsea_df, "q₀" => ByRow(<(0.2))))
    res = fsea(1:row.nfeatures, row.ranks)
    fig = Figure(; size=(800, 800))
    gr = GridLayout(fig[1,1])
    _, ax1, ax2 = plot!(gr, res)

    text!(ax1, 0, 1;
            text = "$(row.timepoint)\n$(row.eeg_feature)\n$(row.geneset)",
            align = (:left, :top),
            offset = (4, -2),
            space = :relative,
            fontsize = 12
    )
    save("data/figures/enrichments/future_$(row.timepoint)_$(row.eeg_feature)_$(replace(row.geneset, " "=>"-")).png", fig)
end

##

tps = ("3m", "6m", "12m")
gdf = groupby(fsea_df, "eeg_feature")
feats = copy(eeg_features)
featidx = Dict(f=> i for (i,f) in enumerate(feats))
gss = unique(fsea_df.geneset)
gsidx = Dict(f=> i for (i,f) in enumerate(gss))

cs = reverse(ColorSchemes.PuOr_11[[1,3,4,6,8,9,11]])
cs[4] = colorant"gray70"

for tp in tps
    @warn "$tp"
    fig = Figure(; size=(1200, 700))
    Legend(fig[2,1:3], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
           ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
            "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
           orientation=:horizontal, tellheight=true, tellwidth=false)
    for (i, feat) in enumerate(filter(contains("amp"), feats))
        @info "$feat"

        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(fig[1,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  title  = replace(feat, "peak_amp_"=>"", "_corrected"=>"c"),
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    Label(fig[0, :], text = "Concurrent Amplitudes; $tp", fontsize = 30)
    save("data/figures/$(tp)_amplitude_fsea_summary.png", fig)
end

for feat in filter(contains("amp"), feats)
    featstr = replace(feat, "peak_"=>"")
    @warn "$featstr"
    fig = Figure(; size=(1200, 700))
    Legend(fig[2,1:3], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
           ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
            "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
           orientation=:horizontal, tellheight=true, tellwidth=false)
    for (i, tp) in enumerate(tps)
        @info "$tp"

        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(fig[1,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  title  = tp,
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    Label(fig[0, :], text = "Concurrent $featstr", fontsize = 30)
    save("data/figures/$(featstr)_fsea_summary.png", fig)
end


for tp in tps
    @warn "$tp"
    fig = Figure(; size=(1200, 700))
    Legend(fig[2,1:3], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
           ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
            "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
           orientation=:horizontal, tellheight=true, tellwidth=false)
    for (i, feat) in enumerate(filter(contains("latency"), feats))
        @info "$feat"

        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(fig[1,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  title  = replace(feat, "peak_latency_"=>"", "_corrected"=>"c"),
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    Label(fig[0, :], text = "Concurrent Latencies; $tp", fontsize = 30)
    save("data/figures/$(tp)_latency_fsea_summary.png", fig)
end

for feat in filter(contains("latency"), feats)
    featstr = replace(feat, "peak_"=>"")
    @warn "$featstr"
    fig = Figure(; size=(1200, 700))
    Legend(fig[2,1:3], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
           ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
            "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
           orientation=:horizontal, tellheight=true, tellwidth=false)
    for (i, tp) in enumerate(tps)
        @info "$tp"

        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(fig[1,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  title  = tp,
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    Label(fig[0, :], text = "Concurrent $featstr", fontsize = 30)
    save("data/figures/$(featstr)_fsea_summary.png", fig)
end

#######

fig = Figure(; size=(1400,1400))
grid = GridLayout(fig[1,1])
Label(grid[-1, 1:3], text = "Concurrent Amplitudes", fontsize=30)

Legend(fig[2,1], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
       ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
        "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
       orientation=:horizontal, tellheight=true, tellwidth=false)
for (j, tp) in enumerate(tps)
    @warn "$tp"
    Label(grid[j,0]; text = tp, fontsize = 20, tellheight=false, tellwidth=true)
    for (i, feat) in enumerate(filter(contains("amp"), feats))
        @info "$feat"

        Label(grid[0,i]; text = replace(feat, r"peak_[a-z]+_"=>"", "_corrected"=>"c"),
              fontsize = 20,
              tellheight=true,
              tellwidth=false)
        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(grid[j,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    save("data/figures/amplitude_fsea_summary.png", fig)
end


fig = Figure(; size=(1400,1400))
grid = GridLayout(fig[1,1])
Label(grid[-1, 1:3], text = "Concurrent Latencies", fontsize=30)

Legend(fig[2,1], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
       ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
        "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
       orientation=:horizontal, tellheight=true, tellwidth=false)
for (j, tp) in enumerate(tps)
    @warn "$tp"
    Label(grid[j,0]; text = tp, fontsize = 20, tellheight=false, tellwidth=true)
    for (i, feat) in enumerate(filter(contains("latency"), feats))
        @info "$feat"

        Label(grid[0,i]; text = replace(feat, r"peak_[a-z]+_"=>"", "_corrected"=>"c"),
              fontsize = 20,
              tellheight=true,
              tellwidth=false)
        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(grid[j,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    save("data/figures/latency_fsea_summary.png", fig)
end


# latency
gdf = groupby(fsea_df, ["eeg_feature", "timepoint"])
lats = filter(f-> contains(f, "latency"), feats)
amps = filter(f-> contains(f, "amp"), feats)
latesmat = zeros(length(gss), 9)
ampesmat = zeros(length(gss), 9)

latqmat = zeros(length(gss), 9)
ampqmat = zeros(length(gss), 9)

for (i, (feat, tp)) in enumerate(Iterators.product(amps, tps))
    df = gdf[(; eeg_feature=feat, timepoint=tp)]
    for row in eachrow(df)
        ampqmat[gsidx[row.geneset], i] = log(row.q₀ + 1e-4) * sign(row.es)
        if row.q₀ < 0.2
            ampesmat[gsidx[row.geneset], i] = row.es
        end
    end
end

    
for (i, (feat, tp)) in enumerate(Iterators.product(lats, tps))
    df = gdf[(; eeg_feature=feat, timepoint=tp)]
    for row in eachrow(df)
        latqmat[gsidx[row.geneset], i] = log(row.q₀ + 1e-4) * sign(row.es)
        if row.q₀ < 0.2
            latesmat[gsidx[row.geneset], i] = row.es
        end
    end
end

ampcl = hclust(pairwise(Euclidean(), ampesmat'); branchorder=:optimal)
atcl = hclust(pairwise(Euclidean(), latesmat'); branchorder=:optimal)
#- 

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, amps), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[ampcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, amps), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, amps), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, ampesmat'[1:3, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, ampesmat'[4:6, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, ampesmat'[7:9, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_amp_heatmap.png", figure)
#-
figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, lats), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[latcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, lats), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, latesmat'[1:3, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, latesmat'[4:6, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, latesmat'[7:9, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_lat_heatmap.png", figure)
#-

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, amps), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[ampcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, amps), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, amps), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, ampqmat'[1:3, ampcl.order]; colormap=:vik, colorrange = (-4.0, 4.0), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, ampqmat'[4:6, ampcl.order]; colormap=:vik, colorrange = (-4.0, 4.0), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, ampqmat'[7:9, ampcl.order]; colormap=:vik, colorrange = (-4.0, 4.0), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="log(qvalue)", )
save("data/figures/concurrent_amp_q_heatmap.png", figure)

#-

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, lats), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[latcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, lats), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, latqmat'[1:3, latcl.order]; colormap=:vik, colorrange = (-4.0, 4.0), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, latqmat'[4:6, latcl.order]; colormap=:vik, colorrange = (-4.0, 4.0), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, latqmat'[7:9, latcl.order]; colormap=:vik, colorrange = (-4.0, 4.0), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="log(qvalue)", )
save("data/figures/concurrent_lat_q_heatmap.png", figure)

# Nodiff
gdf = groupby(fsea_df_nodiff, ["eeg_feature", "timepoint"])
lats = filter(f-> contains(f, "latency"), feats)
amps = filter(f-> contains(f, "amp"), feats)
latesmat = zeros(length(gss), 9)
ampesmat = zeros(length(gss), 9)

latqmat = zeros(length(gss), 9)
ampqmat = zeros(length(gss), 9)

for (i, (feat, tp)) in enumerate(Iterators.product(amps, tps))
    df = gdf[(; eeg_feature=feat, timepoint=tp)]
    for row in eachrow(df)
        ampqmat[gsidx[row.geneset], i] = log(row.q₀ + 1e-4) * sign(row.es)
        if row.q₀ < 0.2
            ampesmat[gsidx[row.geneset], i] = row.es
        end
    end
end

    
for (i, (feat, tp)) in enumerate(Iterators.product(lats, tps))
    df = gdf[(; eeg_feature=feat, timepoint=tp)]
    for row in eachrow(df)
        latqmat[gsidx[row.geneset], i] = log(row.q₀ + 1e-4) * sign(row.es)
        if row.q₀ < 0.2
            latesmat[gsidx[row.geneset], i] = row.es
        end
    end
end

ampcl = hclust(pairwise(Euclidean(), ampesmat'); branchorder=:optimal)
latcl = hclust(pairwise(Euclidean(), latesmat'); branchorder=:optimal)
#- 

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, amps), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[ampcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, amps), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, amps), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, ampesmat'[1:3, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, ampesmat'[4:6, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, ampesmat'[7:9, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_nodiff_amp_heatmap.png", figure)
#-
figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, lats), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[latcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, lats), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, latesmat'[1:3, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, latesmat'[4:6, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, latesmat'[7:9, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_nodiff_lat_heatmap.png", figure)


### Future ####





feats = copy(eeg_features)
featidx = Dict(f=> i for (i,f) in enumerate(feats))
gss = unique(vcat(future6m_fsea_df, future12m_fsea_df).geneset)
gsidx = Dict(f=> i for (i,f) in enumerate(gss))
tps = String.(unique(vcat(future6m_fsea_df, future12m_fsea_df).timepoint))


gdf = groupby(vcat(future6m_fsea_df, future12m_fsea_df), ["eeg_feature", "timepoint"])

fig = Figure(; size=(1400,1400))
grid = GridLayout(fig[1,1])
Label(grid[-1, 1:3], text = "Future Amplitudes", fontsize=30)

Legend(fig[2,1], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
       ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
        "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
       orientation=:horizontal, tellheight=true, tellwidth=false)
for (j, tp) in enumerate(tps)
    @warn "$tp"
    Label(grid[j,0]; text = replace(tp, "_future"=>"->"), fontsize = 20, tellheight=false, tellwidth=true)
    for (i, feat) in enumerate(filter(contains("amp"), feats))
        @info "$feat"

        Label(grid[0,i]; text = replace(feat, r"peak_[a-z]+_"=>"", "_corrected"=>"c"),
              fontsize = 20,
              tellheight=true,
              tellwidth=false)
        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(grid[j,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = gdf[(; eeg_feature=feat, timepoint=tp)]
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    save("data/figures/amplitude_future_fsea_summary.png", fig)
end


fig = Figure(; size=(1400,1400))
grid = GridLayout(fig[1,1])
Label(grid[-1, 1:3], text = "Future Latencies", fontsize=30)

Legend(fig[2,1], [MarkerElement(; marker=:rect, color=cs[i]) for i in [1:3;5:7]],
       ["(-) q < 0.01", "(-) q < 0.1", "(-) q < 0.2", 
        "(+) q < 0.2", "(+) q < 0.1", "(+) q < 0.01"];
       orientation=:horizontal, tellheight=true, tellwidth=false)
for (j, tp) in enumerate(tps)
    @warn "$tp"
    Label(grid[j,0]; text = replace(tp, "_future"=>"->"), fontsize = 20, tellheight=false, tellwidth=true)
    for (i, feat) in enumerate(filter(contains("latency"), feats))
        @info "$feat"

        Label(grid[0,i]; text = replace(feat, r"peak_[a-z]+_"=>"", "_corrected"=>"c"),
              fontsize = 20,
              tellheight=true,
              tellwidth=false)
        df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
        ax = Axis(grid[j,i];
                  ylabel = "geneset",
                  xlabel = "z",
                  yticks = (1:length(gss), gss)
        )
        allymed = median(df.z)
        subdf = gdf[(; eeg_feature=feat, timepoint=tp)]
        lines!(ax, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        hideydecorations!(ax, grid=false, ticks=false, ticklabels = i != 1)
        hidexdecorations!(ax, ticks=false, ticklabels=false)
        for row in eachrow(subdf)
            yidx = findall(u-> replace(u, "UniRef90_"=>"") ∈ na_map[row.geneset], df.feature)
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(yidx)) .+ xpos
            ys = df.z[yidx]
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? cs[cidx] : cs[8 - cidx]
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2), orientation=:horizontal)
            scatter!(ax, ys, xs; color = c)    
            lines!(ax, [ymed, ymed], [xpos - 0.4, xpos + 0.4], color = c)
        end
    end
    save("data/figures/latency_future_fsea_summary.png", fig)
end

lats = filter(f-> contains(f, "latency"), feats)
amps = filter(f-> contains(f, "amp"), feats)
latesmat = zeros(length(gss), 9)
ampesmat = zeros(length(gss), 9)

latqmat = zeros(length(gss), 9)
ampqmat = zeros(length(gss), 9)

for (i, (feat, tp)) in enumerate(Iterators.product(amps, String.(unique(DataFrames.select(gdf, :).timepoint))))
    df = gdf[(; eeg_feature=feat, timepoint=tp)]
    for row in eachrow(df)
        ampqmat[gsidx[row.geneset], i] = log(row.q₀ + 1e-4) * sign(row.es)
        if row.q₀ < 0.2
            ampesmat[gsidx[row.geneset], i] = row.es
        end
    end
end

    
for (i, (feat, tp)) in enumerate(Iterators.product(lats, unique(DataFrames.select(gdf, :).timepoint)))
    df = gdf[(; eeg_feature=feat, timepoint=tp)]
    for row in eachrow(df)
        latqmat[gsidx[row.geneset], i] = log(row.q₀ + 1e-4) * sign(row.es)
        if row.q₀ < 0.2
            latesmat[gsidx[row.geneset], i] = row.es
        end
    end
end

ampcl = hclust(pairwise(Euclidean(), ampesmat'); branchorder=:optimal)
latcl = hclust(pairwise(Euclidean(), latesmat'); branchorder=:optimal)
#- 

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m->6m", xticks = (1:3, amps), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[ampcl.order]))
ax2 = Axis(figure[1,2]; title="3m->12m", xticks = (1:3, amps), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="6m->12m", xticks = (1:3, amps), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, ampesmat'[1:3, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, ampesmat'[4:6, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, ampesmat'[7:9, ampcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/future_amp_heatmap.png", figure)
#-
figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m->6m", xticks = (1:3, lats), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[latcl.order]))
ax2 = Axis(figure[1,2]; title="3m->12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="6m->12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, latesmat'[1:3, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
heatmap!(ax2, latesmat'[4:6, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
hm = heatmap!(ax3, latesmat'[7:9, latcl.order]; colormap=:vik, colorrange = (-0.5, 0.5), highclip=:orangered, lowclip=:indigo)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/future_lat_heatmap.png", figure)

### Bug Contributions ###
using Preferences
using Microbiome

na_unirefs = Set(reduce(union, values(na_map)))

mbotps = Dict(
    "3m"  => concurrent_3m,
    "6m"  => concurrent_6m,
    "12m" => concurrent_12m,
    "3m_future6m" => future_3m6m,
    "3m_future12m"   => future_3m12m,
    "6m_future12m"   => future_6m12m
)
seqs = Set(mapreduce(samplenames, vcat, values(mbotps)))
#


humann_files = mapreduce(vcat, readdir(joinpath(load_preference(VKCComputing, "mgx_analysis_dir"), "humann", "main"); join=true)) do f
    m = match(r"(SEQ\d+)", f)
    isnothing(m) && return DataFrame()
    sample = replace(basename(f), r"(SEQ\d+)_S\d+.+" => s"\1")
    sample ∈ seqs || return DataFrame()
    if contains(basename(f), "genefamilies.tsv") && m[1] ∈ seqs
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

topgenera = mapreduce(vcat, eachrow(vcat(subset(fsea_df, "q₀" => ByRow(<(0.2))),
                                         subset(future6m_fsea_df, "q₀" => ByRow(<(0.2))),
                                         subset(future12m_fsea_df, "q₀" => ByRow(<(0.2)))
            ))) do row
    tp = row.timepoint
    gs = row.geneset
    eeg_feat = row.eeg_feature

    unirefs = na_map[gs]

    df = subset(humann_files, "uniref" => ByRow(u -> u ∈ unirefs))
    df = grouptop(df, 10)
    df = sort(DataFrames.combine(groupby(df, "genus"), "abundance"=>sum => "abundance"), "abundance"; rev=true)

    df.timepoint .= tp
    df.geneset .= gs
    df.eeg_feature .= eeg_feat
    df
end

topspecies = mapreduce(vcat, eachrow(vcat(subset(fsea_df, "q₀" => ByRow(<(0.2))),
                                         subset(future6m_fsea_df, "q₀" => ByRow(<(0.2))),
                                         subset(future12m_fsea_df, "q₀" => ByRow(<(0.2)))
            ))) do row
    tp = row.timepoint
    gs = row.geneset
    eeg_feat = row.eeg_feature

    unirefs = na_map[gs]

    df = subset(humann_files, "uniref" => ByRow(u -> u ∈ unirefs))
    df = grouptop(df, 20; groupcol="species")
    df = sort(DataFrames.combine(groupby(df, "species"), "abundance"=>sum => "abundance"), "abundance"; rev=true)

    df.timepoint .= tp
    df.geneset .= gs
    df.eeg_feature .= eeg_feat
    df
end

CSV.write("data/outputs/fsea_top_genera.csv", topgenera)
CSV.write("data/outputs/fsea_top_species.csv", topspecies)


open("data/outputs/topspecies.txt", "w") do io
    println.(io, filter(!=("other"), unique(topspecies.species)));
end;
open("data/outputs/topgenera.txt", "w") do io
    println.(io, filter(!=("other"), unique(topgenera.genus)))
end;

#-

mbotps_df = Dict(k => begin
                     comm = mbotps[k]
                     df = DataFrame(get(comm))
                     for f in features(comm)
                         df[!, string(f)] = vec(abundances(comm[f,:]))
                     end
                     df
                 end for k in keys(mbotps)
)


for row in eachrow(vcat(subset(fsea_df, "q₀" => ByRow(<(0.2))),
                                         subset(future6m_fsea_df, "q₀" => ByRow(<(0.2))),
                                         subset(future12m_fsea_df, "q₀" => ByRow(<(0.2)))
    ))
    tp = row.timepoint
    gs = row.geneset
    eeg_feat = row.eeg_feature

    unirefs = na_map[gs]

    subdf = let df = subset(humann_files, "uniref" => ByRow(u -> u ∈ unirefs))
        df = leftjoin(mbotps_df[tp], grouptop(df); on="sample")
        subset(df, "genus"=> ByRow(!ismissing))
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
    
#-##############
# Spec cormaps #
################

#### Concurrent ####

gdf = groupby(fsea_df, ["eeg_feature", "timepoint"])

bugmapdir = joinpath("data", "figures", "bugmaps")
isdir(bugmapdir) || mkdir(bugmapdir)

tps = String.(unique(fsea_df.timepoint))
gsp = groupby(topspecies, ["geneset", "eeg_feature"])

for feat in eeg_features, tp in tps 
    @warn feat, tp
    lms = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
    for gs in String.(unique(mapreduce(g-> select(g, "geneset").geneset, union, gdf)))
        sublms = sort(subset(lms, "feature"=> ByRow(f-> replace(f, "UniRef90_"=>"") ∈ na_map[gs])), "z")
        subsp = sort(subset(get(gsp, (; geneset=gs, eeg_feature=feat), DataFrame(timepoint=String[], abundance=Float64[])),
                            "timepoint"=> ByRow(==(tps[3]))), "abundance"; rev=true)
        any(isempty, (sublms, subsp)) && continue
        @info gs

        subhum = Dict()
        for row in eachrow(subset(DataFrames.combine(groupby(humann_files, "feature"),
                                           "abundance"=> sum => "abundance",
                                           "uniref"=> first => "uniref",
                                           "species"=> first=> "species"),
                        "uniref"=> ByRow(u-> u ∈ na_map[gs])
            ))
            subhum[row.species] = get(subhum, row.species, Dict())
            subhum[row.species][row.uniref] = row.abundance
        end
        mat = zeros(size(subsp, 1), size(sublms, 1))

        allsp = Set(filter(!=("other"), subsp.species))

        for (i, sp) in enumerate(subsp.species), (j, uniref) in enumerate(sublms.feature)
            uni = replace(uniref, "UniRef90_"=>"")
            if sp == "other"
                mat[i, j] = mapreduce(+, keys(subhum)) do sp
                    sp ∈ allsp && return 0.0
                    get(subhum[sp], uni, 0.0)
                end
            else
                mat[i, j] = get(subhum[sp], uni, 0.0)
            end
        end

        fig = Figure(; size=(1000,500))
        grid = GridLayout(fig[1,1])
        ax_hm = Axis(grid[1,1]; title = string(feat, "; ", tp, "; ", gs))
        hidedecorations!(ax_hm)
        ax_u = Axis(grid[1,0]; ylabel="$gs genes")
        hideydecorations!(ax_u; label=false)
        hidexdecorations!(ax_u)
        ax_bug = Axis(grid[2,1]; xlabel="species", xticks = (1:size(subsp, 1), subsp.species), xticklabelrotation=π/4)
        hideydecorations!(ax_bug)

        hm = heatmap!(ax_hm, mat)
        uhm = heatmap!(ax_u, reshape(sublms.z, 1, size(sublms, 1)); colormap=:vik)
        bughm = heatmap!(ax_bug, log.(reshape(subsp.abundance, size(subsp, 1), 1)); colormap= :magma)
        colgap!(grid, 10)
        rowgap!(grid, 10)
        colsize!(grid, 0, Fixed(20))
        rowsize!(grid, 2, Fixed(20))

        Colorbar(fig[1,2], hm; label = "RPKM")
        Colorbar(fig[1,3], uhm; label = "z statistic")
        Colorbar(fig[1,4], bughm; label = "log₂(abundance)")

        save(joinpath(bugmapdir, "$(feat)_$(tp)_$(replace(gs, " " => "-")).png"), fig)
    end
end

#### Future ####

gdf = groupby(vcat(future6m_fsea_df, future12m_fsea_df), ["eeg_feature", "timepoint"])

bugmapdir = joinpath("data", "figures", "bugmaps")
isdir(bugmapdir) || mkdir(bugmapdir)

tps = String.(unique(vcat(future6m_fsea_df, future12m_fsea_df).timepoint))
gsp = groupby(topspecies, ["geneset", "eeg_feature"])

for feat in eeg_features, tp in tps 
    @warn feat, tp
    lms = df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
    for gs in String.(unique(mapreduce(g-> select(g, "geneset").geneset, union, gdf)))
        sublms = sort(subset(lms, "feature"=> ByRow(f-> replace(f, "UniRef90_"=>"") ∈ na_map[gs])), "z")
        subsp = sort(subset(get(gsp, (; geneset=gs, eeg_feature=feat), DataFrame(timepoint=String[], abundance=Float64[])),
                            "timepoint"=> ByRow(==(tps[3]))), "abundance"; rev=true)
        any(isempty, (sublms, subsp)) && continue
        @info gs

        subhum = Dict()
        for row in eachrow(subset(DataFrames.combine(groupby(humann_files, "feature"),
                                           "abundance"=> sum => "abundance",
                                           "uniref"=> first => "uniref",
                                           "species"=> first=> "species"),
                        "uniref"=> ByRow(u-> u ∈ na_map[gs])
            ))
            subhum[row.species] = get(subhum, row.species, Dict())
            subhum[row.species][row.uniref] = row.abundance
        end
        mat = zeros(size(subsp, 1), size(sublms, 1))

        allsp = Set(filter(!=("other"), subsp.species))

        for (i, sp) in enumerate(subsp.species), (j, uniref) in enumerate(sublms.feature)
            uni = replace(uniref, "UniRef90_"=>"")
            if sp == "other"
                mat[i, j] = mapreduce(+, keys(subhum)) do sp
                    sp ∈ allsp && return 0.0
                    get(subhum[sp], uni, 0.0)
                end
            else
                mat[i, j] = get(subhum[sp], uni, 0.0)
            end
        end

        fig = Figure(; size=(1000,500))
        grid = GridLayout(fig[1,1])
        ax_hm = Axis(grid[1,1]; title = string(feat, "; ", tp, "; ", gs))
        hidedecorations!(ax_hm)
        ax_u = Axis(grid[1,0]; ylabel="$gs genes")
        hideydecorations!(ax_u; label=false)
        hidexdecorations!(ax_u)
        ax_bug = Axis(grid[2,1]; xlabel="species", xticks = (1:size(subsp, 1), subsp.species), xticklabelrotation=π/4)
        hideydecorations!(ax_bug)

        hm = heatmap!(ax_hm, mat)
        uhm = heatmap!(ax_u, reshape(sublms.z, 1, size(sublms, 1)); colormap=:vik)
        bughm = heatmap!(ax_bug, log.(reshape(subsp.abundance, size(subsp, 1), 1)); colormap= :magma)
        colgap!(grid, 10)
        rowgap!(grid, 10)
        colsize!(grid, 0, Fixed(20))
        rowsize!(grid, 2, Fixed(20))

        Colorbar(fig[1,2], hm; label = "RPKM")
        Colorbar(fig[1,3], uhm; label = "z statistic")
        Colorbar(fig[1,4], bughm; label = "log₂(abundance)")

        save(joinpath(bugmapdir, "$(feat)_$(tp)_$(replace(gs, " " => "-")).png"), fig)
    end
end

#-
all_3m = commjoin(concurrent_3m, future_3m12m[:, setdiff(samplenames(future_3m12m), samplenames(concurrent_3m))],
                                 future_3m6m[:, setdiff(samplenames(future_3m6m), union(samplenames(concurrent_3m), samplenames(future_3m12m)))])
all_6m = commjoin(concurrent_6m, future_6m12m[:, setdiff(samplenames(future_6m12m), samplenames(concurrent_6m))])
all_12m = copy(concurrent_12m)

#- 

ages = mapreduce(comm-> get(comm, :age), vcat, (all_3m, all_6m, all_12m))
ntop = 8

tp = "3m_future12m"
gs = "GABA synthesis"
feat = "peak_amp_N2_corrected"

sps = subset(topspecies, "timepoint"=> ByRow(==(tp)),
                        "geneset"=> ByRow(==(gs)),
                        "eeg_feature"=> ByRow(==(feat)),
                        "species" => ByRow(!=("other"))
)

colors = [sp=> c for (sp, c) in zip(sps.species[1:ntop], cgrad(:managua, ntop; categorical=true))]

df = mapreduce(vcat, sps.species[1:ntop]) do sp
    comms = filter(c-> haskey(c.fidx, sp), (all_3m, all_6m, all_12m))
    abunds = mapreduce(comm-> vec(collect(abundances(comm[sp, :]))), vcat, comms)
    ixs = findall(>(0), abunds)
    isempty(ixs) && return DataFrame(age=Float64[], abundance=Float64[], species=String[])

    DataFrame(age = ages[ixs], abundance = abunds[ixs], species = fill(sp, length(ixs)))
end

xyz = data(df) * mapping(:age, :abundance, color=:species, layout=:species)
layers = smooth() + visual(Scatter)
fg = draw(layers * xyz, palettes=(; color=colors), figure=(; size=(1500, 1500)))


#text!(12, sum(abundances(all_12m[sp,:])); text = [rich(replace(sp, r"s__([A-Z])[a-z]+_"=>s"\1. "); font=:italic, color=:dodgerblue)])  
save("data/figures/bugs/top_$(tp)_$(gs)_$(feat)_line.png", fg)
