using FeatureSetEnrichments
using VKCComputing
using EEGMicrobiome
using DataFrames
using Distances
using CSV
using HypothesisTests
using MultipleTesting
using Clustering
using GLM
using Distributions

# load data
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
transform!(future_3m12m_func, AsTable(["age", "eeg_age"])=> ByRow(nt-> nt.eeg_age - nt.age)=> "age_diff")
load_functional_profiles!(future_3m12m_func);

future_6m12m_func = DataFrame(get(future_6m12m))
transform!(future_6m12m_func, AsTable(["age", "eeg_age"])=> ByRow(nt-> nt.eeg_age - nt.age)=> "age_diff")
load_functional_profiles!(future_6m12m_func);

##


eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency"=>"amp")]

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()
na_map_full = FeatureSetEnrichments.get_neuroactive_unirefs(; consolidate=false)

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(concurrent_3m_func, "./data/outputs/lms/$(feature)_nodiff_3m_lms.csv",
                         feature, names(concurrent_3m_func, r"^UniRef")
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(concurrent_6m_func, "./data/outputs/lms/$(feature)_nodiff_6m_lms.csv",
                         feature, names(concurrent_6m_func, r"^UniRef")
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(concurrent_12m_func, "./data/outputs/lms/$(feature)_nodiff_12m_lms.csv",
                         feature, names(concurrent_12m_func, r"^UniRef")
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(concurrent_3m_func, "./data/outputs/lms/$(feature)_3m_lms.csv",
                         feature, names(concurrent_3m_func, r"^UniRef")
        additional_cols = [:age_diff],
        formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(concurrent_6m_func, "./data/outputs/lms/$(feature)_6m_lms.csv",
                         feature, names(concurrent_6m_func, r"^UniRef")
        additional_cols = [:age_diff],
        formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(concurrent_12m_func, "./data/outputs/lms/$(feature)_12m_lms.csv",
                         feature, names(concurrent_12m_func, r"^UniRef")
        additional_cols = [:age_diff],
        formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(future_3m6m_func, "./data/outputs/lms/$(feature)_3m_future6m_lms.csv", feature, names(future_3m6m_func, r"^UniRef");
        additional_cols = [:age_diff],
        formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(future_3m12m_func, "./data/outputs/lms/$(feature)_3m_future12m_lms.csv", feature, names(future_3m12m_func, r"^UniRef");
        additional_cols = [:age_diff],
        formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
    )
end

for feature in eeg_features
    @info feature
    EEGMicrobiome.runlms(future_6m12m_func, "./data/outputs/lms/$(feature)_6m_future12m_lms.csv", feature, names(future_6m12m_func, r"^UniRef");
        additional_cols = [:age_diff],
        formula = term(:func) ~ term(feature) + term(:age) + term(:age_diff) + term(:n_segments)
    )
end

## Run FSEA

tps = ("3m", "6m", "12m")

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
# fsea_df = CSV.read("data/outputs/fsea/true_ages_fsea.csv", DataFrame)

fsea_df_nodiff = let fsea_df_nodiff = DataFrame()
    for tp in tps, feature in eeg_features
        filepath = "./data/outputs/lms/$(feature)_nodiff_$(tp)_lms.csv"
        lms = CSV.read(filepath, DataFrame)

        for (key, unirefs) in pairs(na_map)
            s = Set("UniRef90_$uniref" for uniref in unirefs)
            lms[!, key] = lms.feature .∈ Ref(s)
        end

        for (key, unirefs) in pairs(na_map)
            idx = findall(lms[!, key])
            length(idx) > 5 || continue
            result = fsea(FeatureSetEnrichments.Permutation(5000), lms.z, idx)
            push!(fsea_df_nodiff, (; timepoint=tp,
                              eeg_feature = feature,
                              geneset     = key,
                              pvalue      = pvalue(result),
                              es          = enrichment_score(result),
                              ranks       = idx,
                              nfeatures   = size(lms, 1)))
        end
    end
    fsea_df_nodiff
end

transform!(fsea_df_nodiff, "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "q₀")
transform!(groupby(fsea_df_nodiff, "timepoint"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qₜ")
transform!(groupby(fsea_df_nodiff, "geneset"), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg()))=> "qᵧ")
sort!(fsea_df_nodiff, "q₀")
CSV.write("data/outputs/fsea/concurrent_nodiff_consolidated_fsea.csv", fsea_df_nodiff)
# fsea_df = CSV.read("data/outputs/fsea/true_ages_fsea.csv", DataFrame)
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
            length(idx) > 5 || continue
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
CSV.write("data/outputs/fsea/future6m_consolidated_fsea.csv", future6m_fsea_df)
# future6m_fsea_df = CSV.read("data/outputs/fsea/future6m_true_ages_fsea.csv", DataFrame)


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

gdf = groupby(fsea_df, "eeg_feature")
feats = copy(eeg_features)
featidx = Dict(f=> i for (i,f) in enumerate(feats))
gss = unique(fsea_df.geneset)
gsidx = Dict(f=> i for (i,f) in enumerate(gss))

for feat in feats
    for tp in tps
        fig = Figure(; size=(700, 500))
        ax = Axis(fig[1,1], xlabel="geneset", ylabel = "rank", title=feat, xticks=(1:length(gss), gss), xticklabelrotation = pi/4)
        for row in eachrow(subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp))))
            xpos = gsidx[row.geneset]
            xs = rand(Normal(0.0, 0.05), length(row.ranks)) .+ xpos
            med = median(1:row.nfeatures)
            ys = (row.ranks .- med) ./ med
            ymed = median(ys)
            c = row.q₀ > 0.2 ? :gray : ymed < 0 ? :blue : :red
            violin!(ax, fill(xpos, length(ys)), ys; color=(c, 0.2))
            scatter!(ax, xs, ys; color = c)    
            ylims!(ax, (-1.5, 1.5))
            #lines!(ax, [xpos - 0.5, xpos + 0.5], [ymed, ymed], color = c)
        end
        save("data/figures/$(feat)_$(tp)_fsea_summary.png", fig)
    end
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
latcl = hclust(pairwise(Euclidean(), latesmat'); branchorder=:optimal)
#- 

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, amps), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[ampcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, amps), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, amps), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, ampesmat'[1:3, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
heatmap!(ax2, ampesmat'[4:6, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
hm = heatmap!(ax3, ampesmat'[7:9, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_amp_heatmap.png", figure)
#-
figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, lats), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[latcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, lats), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, latesmat'[1:3, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
heatmap!(ax2, latesmat'[4:6, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
hm = heatmap!(ax3, latesmat'[7:9, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_lat_heatmap.png", figure)
#-

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, amps), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[ampcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, amps), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, amps), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, ampqmat'[1:3, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-4.0, 4.0), highclip=:yellow, lowclip=:gray10)
heatmap!(ax2, ampqmat'[4:6, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-4.0, 4.0), highclip=:yellow, lowclip=:gray10)
hm = heatmap!(ax3, ampqmat'[7:9, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-4.0, 4.0), highclip=:yellow, lowclip=:gray10)
Colorbar(figure[1,4], hm; label="log(qvalue)", )
save("data/figures/concurrent_amp_q_heatmap.png", figure)

#-

figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, lats), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[latcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, lats), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, latqmat'[1:3, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-4.0, 4.0), highclip=:yellow, lowclip=:gray10)
heatmap!(ax2, latqmat'[4:6, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-4.0, 4.0), highclip=:yellow, lowclip=:gray10)
hm = heatmap!(ax3, latqmat'[7:9, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-4.0, 4.0), highclip=:yellow, lowclip=:gray10)
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

heatmap!(ax1, ampesmat'[1:3, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
heatmap!(ax2, ampesmat'[4:6, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
hm = heatmap!(ax3, ampesmat'[7:9, ampcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_nodiff_amp_heatmap.png", figure)
#-
figure = Figure(;size = (1200, 800))
ax1 = Axis(figure[1,1]; title="3m", xticks = (1:3, lats), xticklabelrotation=pi/4, yticks = (1:length(gss), gss[latcl.order]))
ax2 = Axis(figure[1,2]; title="6m", xticks = (1:3, lats), xticklabelrotation=pi/4)
ax3 = Axis(figure[1,3]; title="12m", xticks = (1:3, lats), xticklabelrotation=pi/4)
hideydecorations!.([ax2,ax3])

heatmap!(ax1, latesmat'[1:3, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
heatmap!(ax2, latesmat'[4:6, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
hm = heatmap!(ax3, latesmat'[7:9, latcl.order]; colormap=Reverse(:RdBu), colorrange = (-0.5, 0.5), highclip=:yellow, lowclip=:gray10)
Colorbar(figure[1,4], hm; label="E.S.", )
save("data/figures/concurrent_nodiff_lat_heatmap.png", figure)

