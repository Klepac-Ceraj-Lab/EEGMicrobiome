using EEGMicrobiome
using ThreadsX
using Chain
using VKCComputing
using Distributions
using CSV
using DataFrames
using CairoMakie
using Microbiome
using BiobakeryUtils
using CairoMakie
using SparseArrays
using MLJ
using MLJDecisionTreeInterface
using GLM
using StableRNGs

concurrent_3m = load_cohort("concurrent_3m")
concurrent_6m = load_cohort("concurrent_6m")
concurrent_12m = load_cohort("concurrent_12m")
future_3m6m = load_cohort("future_3m6m")
future_3m12m = load_cohort("future_3m12m")
future_6m12m = load_cohort("future_6m12m")

topsp = readlines("data/outputs/topspecies.txt")
topgn = readlines("data/outputs/topgenera.txt")

eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency"=>"amp")]

mbotps = Dict(
    "3m"  => concurrent_3m,
    "6m"  => concurrent_6m,
    "12m" => concurrent_12m,
    "3m_future6m" => future_3m6m,
    "3m_future12m"   => future_3m12m,
    "6m_future12m"   => future_6m12m
)


mbotps_df = Dict(k => begin
                     comm = mbotps[k]
                     df = DataFrame(get(comm))
		     transform!(df, AsTable(["age", "age_diff"]) => ByRow(nt-> nt.age + nt.age_diff) => "eeg_age")
                     for f in features(comm)
                         df[!, string(f)] = vec(abundances(comm[f,:]))
                     end
                     df
                 end for k in keys(mbotps)
)

#-

rng_isaiah = StableRNG(2019)

tree = let
    Tree = @load RandomForestRegressor pkg=DecisionTree verbosity=0
    Tree(; n_trees = 500, rng = rng_isaiah)
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

rng_samara = StableRNG(2023)
concurrent_summary_stats = DataFrame()

let modelin = vcat(mbotps_df["3m"], mbotps_df["6m"], mbotps_df["12m"]; cols=:union)
    for col in names(modelin, r"s__")
	modelin[!, col] = coalesce.(modelin[!, col], 0.0)
    end
    for eeg_feat in eeg_features
	@info "Testing $eeg_feat"
	pfolds = partition(unique(modelin.subject), 0.25, 0.25, 0.25; rng=rng_samara)
	trainfolds = map(ps -> findall(s-> s ∉ ps, modelin.subject), pfolds)
	testfolds = map(ps -> findall(s-> s ∈ ps, modelin.subject), pfolds)

	predictions = select(modelin, Not(r"s__"))
	
	for (i, (train, test)) in enumerate(zip(trainfolds, testfolds))
	    @info "Fold $i"
	    age_only_tuned = machine(self_tuning_tree, select(modelin, "eeg_age", "n_segments"), modelin[!, eeg_feat]);
	    plus_bugs_tuned = machine(self_tuning_tree,
				      select(modelin, "eeg_age", "n_segments", intersect(topsp, names(modelin))...),
				      modelin[!, eeg_feat]);
	    bugs_only_tuned = machine(self_tuning_tree, select(modelin, "n_segments", intersect(topsp, names(modelin))...),
				      modelin[!, eeg_feat]);

	    MLJ.fit!(age_only_tuned; rows=train)
	    MLJ.fit!(plus_bugs_tuned; rows=train)
	    MLJ.fit!(bugs_only_tuned; rows=train)

	    age_only_correction = lm(@formula(y ~ yhat + 1),
				     DataFrame("y" => modelin[train, eeg_feat],
					       "yhat" => MLJ.predict(age_only_tuned; rows=train)
					       )
				     )
	    plus_bugs_correction = lm(@formula(y ~ yhat + 1),
				     DataFrame("y" => modelin[train, eeg_feat],
					       "yhat" => MLJ.predict(plus_bugs_tuned; rows=train)
					       )
				     )
	    bugs_only_correction = lm(@formula(y ~ yhat + 1),
				     DataFrame("y" => modelin[train, eeg_feat],
					       "yhat" => MLJ.predict(bugs_only_tuned; rows=train)
					       )
				     )

	    age_only_predictions = MLJ.predict(age_only_tuned; rows=test)
	    plus_bugs_predictions = MLJ.predict(plus_bugs_tuned; rows=test)
	    bugs_only_predictions =  MLJ.predict(bugs_only_tuned; rows=test)
	    #age_only_predictions = GLM.predict(age_only_correction, DataFrame("yhat" => MLJ.predict(age_only_tuned; rows=test)))
	    #plus_bugs_predictions = GLM.predict(plus_bugs_correction, DataFrame("yhat" => MLJ.predict(plus_bugs_tuned; rows=test)))
	    #bugs_only_predictions = GLM.predict(bugs_only_correction, DataFrame("yhat" => MLJ.predict(bugs_only_tuned; rows=test)))

	    age_only_dfinsert  = Vector{Union{Missing, Float64}}(missing, size(modelin, 1))
	    plus_bugs_dfinsert = Vector{Union{Missing, Float64}}(missing, size(modelin, 1))
	    bugs_only_dfinsert = Vector{Union{Missing, Float64}}(missing, size(modelin, 1))
 
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
	    append!(concurrent_summary_stats, DataFrame(
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

	save("data/figures/rf_models/concurrent_$(eeg_feat)_all_ages.png", fig)
    end
end

CSV.write("data/outputs/rf_concurrent_fold_stats.csv", concurrent_summary_stats)
CSV.write("data/outputs/rf_concurren_summary_stats.csv",
