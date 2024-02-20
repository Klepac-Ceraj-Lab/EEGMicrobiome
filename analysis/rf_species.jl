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

eeg = load_eeg()
eeg.peak_latency_P1_corrected = eeg.peak_latency_P1 .- eeg.peak_latency_N1
eeg.peak_latency_N2_corrected = eeg.peak_latency_N2 .- eeg.peak_latency_P1
eeg.peak_amp_P1_corrected = eeg.peak_amp_P1 .- eeg.peak_amp_N1
eeg.peak_amp_N2_corrected = eeg.peak_amp_N2 .- eeg.peak_amp_P1
rename!(eeg, "age" => "eeg_age")
mbo = load_microbiome(eeg.subject)
subset!(mbo, :visit=> ByRow(!ismissing))
transform!(mbo, "visit" => ByRow(v -> ismissing(v) ? v : replace(v, "mo" => "m")) => "timepoint")
load_taxonomic_profiles!(mbo)

topspec = readlines("data/outputs/topspecies.txt")

modelin = @chain eeg begin
    leftjoin(mbo; on=["subject", "timepoint"])
    dropmissing!
    subset!("seqprep" => ByRow(!ismissing))
end
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

summary_stats = DataFrame()

for eeg_feat in names(modelin, r"peak_")
    @info "Testing $eeg_feat"
    pfolds = partition(unique(modelin.subject), 0.25, 0.25, 0.25)
    trainfolds = map(ps -> findall(s-> s ∉ ps, modelin.subject), pfolds)
    testfolds = map(ps -> findall(s-> s ∈ ps, modelin.subject), pfolds)

    predictions = copy(modelin)
    
    for (i, (train, test)) in enumerate(zip(trainfolds, testfolds))
        @info "Fold $i"
        age_only_pred = Union{Missing, Float64}[missing for _ in 1:size(modelin,1)]
        age_only_tuned = machine(self_tuning_tree, select(modelin, "eeg_age", "n_segments"), modelin[!, eeg_feat]);
        plus_bugs_tuned = machine(self_tuning_tree, select(modelin, "eeg_age", "n_segments", topspec...), modelin[!, eeg_feat]);
        bugs_only_tuned = machine(self_tuning_tree, select(modelin, "n_segments", topspec...), modelin[!, eeg_feat]);

        MLJ.fit!(age_only_tuned; rows=train)
        MLJ.fit!(plus_bugs_tuned; rows=train)
        MLJ.fit!(bugs_only_tuned; rows=train)

        age_only_predictions = MLJ.predict(age_only_tuned; rows=test)
        plus_bugs_predictions = MLJ.predict(plus_bugs_tuned; rows=test)
        bugs_only_predictions = MLJ.predict(bugs_only_tuned; rows=test)

        age_only_dfinsert  = Union{Missing, Float64}[missing for _ in 1:size(modelin, 1)]
        plus_bugs_dfinsert = Union{Missing, Float64}[missing for _ in 1:size(modelin, 1)]
        bugs_only_dfinsert = Union{Missing, Float64}[missing for _ in 1:size(modelin, 1)]

        age_only_dfinsert[test]  .= age_only_predictions 
        plus_bugs_dfinsert[test] .= plus_bugs_predictions
        bugs_only_dfinsert[test] .= bugs_only_predictions
        
        predictions[!, "age_only_fold$i"] = age_only_dfinsert
        predictions[!, "plus_bugs_fold$i"] = plus_bugs_dfinsert
        predictions[!, "bugs_only_fold$i"] = bugs_only_dfinsert
         
        age_only_err = mape(age_only_predictions, modelin[test, eeg_feat])
        age_only_r² = cor(age_only_predictions, modelin[test, eeg_feat])^2

        plus_bugs_err = mape(plus_bugs_predictions, modelin[test, eeg_feat])
        plus_bugs_r² = cor(plus_bugs_predictions, modelin[test, eeg_feat])^2

        bugs_only_err = mape(bugs_only_predictions, modelin[test, eeg_feat])
        bugs_only_r² = cor(bugs_only_predictions, modelin[test, eeg_feat])^2
        append!(summary_stats, DataFrame(
             feature = fill(eeg_feat, 3),
             model = ["age_only", "bugs_only", "age_plus_bugs"],
             r² = [age_only_r², bugs_only_r², plus_bugs_r²],
             mape = [age_only_err, bugs_only_err, plus_bugs_err],
             fold = fill(i, 3),
        ))
    end
    
    transform!(predictions, 
        AsTable(r"age_only_fold")=> ByRow(nt-> mean(skipmissing(values(nt))))=> "age_only_mean_pred",
        AsTable(r"plus_bugs_fold")=> ByRow(nt-> mean(skipmissing(values(nt))))=> "plus_bugs_mean_pred"
    )    
    age_only_predictions = predictions.age_only_mean_pred
    plus_bugs_predictions = predictions.plus_bugs_mean_pred

    age_only_err = mape(age_only_predictions, modelin[!, eeg_feat])
    age_only_r² = cor(age_only_predictions, modelin[!, eeg_feat])^2

    plus_bugs_err = mape(plus_bugs_predictions, modelin[!, eeg_feat])
    plus_bugs_r² = cor(plus_bugs_predictions, modelin[!, eeg_feat])^2

    fig = Figure(; size=(400, 600))
    ax1 = Axis(fig[1,1];
        ylabel="Predicted $(replace(eeg_feat, "peak_" => ""))",
        xlabel="Observed $(replace(eeg_feat, "peak_" => ""))",
        title="Age only model"
    )
    ax2 = Axis(fig[2,1];
        ylabel="Predicted $(replace(eeg_feat, "peak_" => ""))",
        xlabel="Observed $(replace(eeg_feat, "peak_" => ""))",
        title="Age + bugs model"
    )

    sc1 = scatter!(ax1, modelin[!, eeg_feat], age_only_predictions;
        color = modelin.eeg_age,
    )

    sc2 = scatter!(ax2,modelin[!, eeg_feat], plus_bugs_predictions;
        color = modelin.eeg_age,
    )

    Colorbar(fig[1:2,2], sc1; label= "age (months)")

    age_only_mn = minimum([modelin[!, eeg_feat]; age_only_predictions])
    age_only_mx = maximum([modelin[!, eeg_feat]; age_only_predictions])
    plus_bugs_mn = minimum([modelin[!, eeg_feat]; plus_bugs_predictions])
    plus_bugs_mx = maximum([modelin[!, eeg_feat]; plus_bugs_predictions])

    lines!(ax1, [age_only_mn, age_only_mx], [age_only_mn, age_only_mx]; color=:gray60, linestyle = :dash)
    lines!(ax2, [plus_bugs_mn, plus_bugs_mx], [plus_bugs_mn, plus_bugs_mx]; color=:gray60, linestyle = :dash)



    text!(ax1, 0, 1;
        text = "R² = $(round(age_only_r², digits=3))\nMAPE = $(round(age_only_err, digits=3))",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 12
    )

    text!(ax2, 0, 1;
        text = "R² = $(round(plus_bugs_r², digits=3))\nMAPE = $(round(plus_bugs_err, digits=3))",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 12
    )

    save("data/figures/rf_models/$(eeg_feat)_all_ages.png", fig)
end

transform!(summary_stats, "model"=> ByRow(m-> begin
    bugs = contains(m, "bugs") ? "+" : "-"
    age = contains(m, "age") ? "+" : "-"

    return (; bugs, age)
end) => ["bugs", "age"])

CSV.write("data/outputs/rf_summary_stats.csv", summary_stats)
