using CSV: partitions
using Statistics
using VKCComputing
using EEGMicrobiome
using FeatureSetEnrichments
using DataFrames
using CSV
using Chain
using CairoMakie
using Preferences
using ColorSchemes
using MLJ
using MLJDecisionTreeInterface
using StableRNGs

tps = ("v1", "v2", "v3")
ftps = ("v1v2", "v1v3", "v2v3")

mdata = load_cohorts()

v1 = get_cohort(mdata, "v1")
v2 = get_cohort(mdata, "v2")
v3 = get_cohort(mdata, "v3")
v1v2 = get_cohort(mdata, "v1v2")
v1v3 = get_cohort(mdata, "v1v3")
v2v3 = get_cohort(mdata, "v2v3")



##

eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency" => "amp")]

na_ko = FeatureSetEnrichments.get_neuroactive_kos()
na_map = FeatureSetEnrichments.get_neuroactive_unirefs()

allkos = Set(union(values(na_ko)...))
allseqs = Set(skipmissing(mdata.seqprep))


kos = let files = filter(f-> contains(basename(f),"kos"), readdir("/grace/sequencing/processed/mgx/humann/rename/", join=true))
    mapreduce(vcat, files) do file
        m = match(r"^([A-Za-z0-9]+)_(S\d+)?_kos_rename\.tsv", basename(file))
        if isnothing(m)
            @error file
            return DataFrame()
        end
        (seq, s_well) = m.captures
        seq in allseqs || return DataFrame()
        @info "loeading $seq"
        df = CSV.read(file, DataFrame; delim='\t', header=["gene", "abundance"], skipto=2)
        df.seqprep .= String(seq)
        return df
    end
end

transform!(kos, "gene"=> ByRow(g-> begin
    spl = split(g, '|')
    (gene, species) = length(spl) > 1 ? spl : (only(spl), missing)
    ko = match(r"K\d+", gene)
    ko = isnothing(ko) ? missing : ko.match
    (; gene, species, ko)
end)=> ["gene", "species", "ko"])

transform!(groupby(kos, "seqid"), ["abundance", "species"]=> ((abundance, species) -> begin
    unstrat = findall(ismissing, species)
    total = sum(abundance[unstrat])
    abundance ./ total .* 100
end) => "relab")

neuro = subset(kos, "ko" => ByRow(ko-> !ismissing(ko) && ko in allkos), "species" => ByRow(ismissing))

neuromat = let mat = unstack(neuro, "seqprep", "ko", "abundance")
    for col in names(mat, r"K\d+")
        mat[!, col] = coalesce.(mat[!, col], 0.)
    end
    mat
end

leftjoin!(mdata, neuromat; on="seqprep", matchmissing=:notequal)



#-


rng = StableRNG(1984)

tree = let
    Tree = @load RandomForestRegressor pkg=DecisionTree verbosity=0
    Tree(; n_trees = 500, rng)
end

r_min_samples_split = range(tree, :min_samples_split; lower=2, upper=10) 
r_max_depth         = range(tree, :max_depth; lower=1, upper=6)
r_min_samples_leaf  = range(tree, :min_samples_leaf; lower=1, upper=5)

self_tuning_tree = TunedModel(
    model = tree,
    tuning = Grid(; resolution=5),
    range = [r_min_samples_split, r_max_depth, r_min_samples_leaf],
    measure = mape
)

#-

modelin = subset(mdata, AsTable(["cohort_v1", "cohort_v2", "cohort_v3"]) => ByRow(any), "K00010"=>ByRow(!ismissing))

for col in names(modelin, r"latency")
    modelin[!, col] = Float64.(modelin[!, col])
end

summary_stats = DataFrame()

for eeg_feat in names(modelin, r"peak_")
    @info "Testing $eeg_feat"
    pfolds = partition(unique(modelin.subject_id), 0.25, 0.25, 0.25)
    trainfolds = map(ps -> findall(s-> s ∉ ps, modelin.subject_id), pfolds)
    testfolds = map(ps -> findall(s-> s ∈ ps, modelin.subject_id), pfolds)

    predictions = copy(modelin)
    
    for (i, (train, test)) in enumerate(zip(trainfolds, testfolds))
        @info "Fold $i"
        age_only_pred = Union{Missing, Float64}[missing for _ in 1:size(modelin,1)]
        age_only_tuned = machine(self_tuning_tree, select(modelin, "eeg_age", "n_trials"), modelin[!, eeg_feat]);
        plus_na_tuned = machine(self_tuning_tree, select(modelin, "eeg_age", "n_trials", Cols(r"K\d+")), modelin[!, eeg_feat]);
        na_only_tuned = machine(self_tuning_tree, select(modelin, "n_trials", Cols(r"K\d+")), modelin[!, eeg_feat]);

        MLJ.fit!(age_only_tuned; rows=train)
        MLJ.fit!(plus_na_tuned; rows=train)
        MLJ.fit!(na_only_tuned; rows=train)

        age_only_predictions = MLJ.predict(age_only_tuned; rows=test)
        plus_na_predictions = MLJ.predict(plus_na_tuned; rows=test)
        na_only_predictions = MLJ.predict(na_only_tuned; rows=test)

        age_only_dfinsert  = Union{Missing, Float64}[missing for _ in 1:size(modelin, 1)]
        plus_na_dfinsert = Union{Missing, Float64}[missing for _ in 1:size(modelin, 1)]
        na_only_dfinsert = Union{Missing, Float64}[missing for _ in 1:size(modelin, 1)]

        age_only_dfinsert[test]  .= age_only_predictions 
        plus_na_dfinsert[test] .= plus_na_predictions
        na_only_dfinsert[test] .= na_only_predictions
        
        predictions[!, "age_only_fold$i"] = age_only_dfinsert
        predictions[!, "plus_na_fold$i"] = plus_na_dfinsert
        predictions[!, "na_only_fold$i"] = na_only_dfinsert
         
        age_only_err = mape(age_only_predictions, modelin[test, eeg_feat])
        age_only_r² = cor(age_only_predictions, modelin[test, eeg_feat])^2

        plus_na_err = mape(plus_na_predictions, modelin[test, eeg_feat])
        plus_na_r² = cor(plus_na_predictions, modelin[test, eeg_feat])^2

        na_only_err = mape(na_only_predictions, modelin[test, eeg_feat])
        na_only_r² = cor(na_only_predictions, modelin[test, eeg_feat])^2
        append!(summary_stats, DataFrame(
             feature = fill(eeg_feat, 3),
             model = ["age_only", "na_only", "age_plus_na"],
             r² = [age_only_r², na_only_r², plus_na_r²],
             mape = [age_only_err, na_only_err, plus_na_err],
             fold = fill(i, 3),
        ))
    end
    
    transform!(predictions, 
        AsTable(r"age_only_fold")=> ByRow(nt-> mean(skipmissing(values(nt))))=> "age_only_mean_pred",
        AsTable(r"plus_na_fold")=> ByRow(nt-> mean(skipmissing(values(nt))))=> "plus_na_mean_pred"
    )

    age_only_predictions = predictions.age_only_mean_pred
    plus_na_predictions = predictions.plus_na_mean_pred

    age_only_err = mape(age_only_predictions, modelin[!, eeg_feat])
    age_only_r² = cor(age_only_predictions, modelin[!, eeg_feat])^2

    plus_na_err = mape(plus_na_predictions, modelin[!, eeg_feat])
    plus_na_r² = cor(plus_na_predictions, modelin[!, eeg_feat])^2

    fig = Figure(; size=(400, 600))
    ax1 = Axis(fig[1,1];
        ylabel="Predicted $(replace(eeg_feat, "peak_" => ""))",
        xlabel="Observed $(replace(eeg_feat, "peak_" => ""))",
        title="Age only model"
    )
    ax2 = Axis(fig[2,1];
        ylabel="Predicted $(replace(eeg_feat, "peak_" => ""))",
        xlabel="Observed $(replace(eeg_feat, "peak_" => ""))",
        title="Age + na model"
    )

    sc1 = scatter!(ax1, modelin[!, eeg_feat], age_only_predictions;
                   color = Float64.(modelin.eeg_age),
    )

    sc2 = scatter!(ax2,modelin[!, eeg_feat], plus_na_predictions;
        color = Float64.(modelin.eeg_age),
    )

    Colorbar(fig[1:2,2], sc1; label= "age (months)")

    age_only_mn = minimum([modelin[!, eeg_feat]; age_only_predictions])
    age_only_mx = maximum([modelin[!, eeg_feat]; age_only_predictions])
    plus_na_mn = minimum([modelin[!, eeg_feat]; plus_na_predictions])
    plus_na_mx = maximum([modelin[!, eeg_feat]; plus_na_predictions])

    lines!(ax1, [age_only_mn, age_only_mx], [age_only_mn, age_only_mx]; color=:gray60, linestyle = :dash)
    lines!(ax2, [plus_na_mn, plus_na_mx], [plus_na_mn, plus_na_mx]; color=:gray60, linestyle = :dash)



    text!(ax1, 0, 1;
        text = "R² = $(round(age_only_r², digits=5))\nMAPE = $(round(age_only_err, digits=3))",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 12
    )

    text!(ax2, 0, 1;
        text = "R² = $(round(plus_na_r², digits=5))\nMAPE = $(round(plus_na_err, digits=3))",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 12
    )

    save("data/figures/rf_models/$(eeg_feat)_neuroactive_all_ages.png", fig)
end

transform!(summary_stats, "model"=> ByRow(m-> begin
    bugs = contains(m, "na") ? "+" : "-"
    age = contains(m, "age") ? "+" : "-"

    return (; bugs, age)
end) => ["na", "age"])

CSV.write("data/outputs/rf_na_summary_stats.csv", summary_stats)
