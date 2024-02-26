using AlgebraOfGraphics: categorical
# # EEG and the Microbiome - figures
# 
# This notebook contains code and explanations
# for plotting the figures (main and supplemental)
# in the manuscript *{Title TBD}*.
#
# ## Setup and code loading
#
# First, we'll load the packages we need.
# The first, `EEGMicrobiome` is the helper code
# written in this repository
# (which can be found in the `src/` directory).
# Other dependecies are encoded in the `Project.toml`
# file at the root of this repo.
#
# Before running this code, you may need to run
# `] instantiate` from the Pkg REPL.

using EEGMicrobiome
using VKCComputing
using FeatureSetEnrichments
using Preferences
using ThreadsX
using Chain
using XLSX
using FileIO
using CSV
sing DataFrames
using Distributions
using Microbiome
using BiobakeryUtils
using CairoMakie
using AlgebraOfGraphics
using SparseArrays
using Clustering
using Distances

# Next, we'll load the data.
# The `load_cohort()` fucntion is written in `src/data_loading.jl`.

tps = ("3m", "6m", "12m")
ftps = ("3m_future6m", "3m_future12m", "6m_future12m")

concurrent_3m = load_cohort("concurrent_3m")
concurrent_6m = load_cohort("concurrent_6m")
concurrent_12m = load_cohort("concurrent_12m")
future_3m6m = load_cohort("future_3m6m")
future_3m12m = load_cohort("future_3m12m")
future_6m12m = load_cohort("future_6m12m")

##

eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency"=>"amp")]

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()


eegmbo = let 
	eeg = load_eeg()
	rename!(eeg, "age"=> "eeg_age")
	mbo = load_microbiome(eeg.subject)
	transform!(mbo, "visit"=>ByRow(v-> replace(v, "mo"=>"m"))=> "visit")

	@chain eeg begin
		select("subject", "timepoint"=>"visit", Cols(:))
		unique!(["subject", "timepoint"])
		outerjoin(mbo; on = ["subject", "visit"])
		groupby("subject")
		subset(AsTable(["visit", "seqprep", "eeg_age"])=> (nt-> begin
			# skip subjects without at least 1 eeg and at least 1 seqprep
			(all(ismissing, nt.seqprep) || all(ismissing(nt.eeg_age))) && return false
			# keep subjects with at least 1 concurrent eeg/seqprep
			any(.!ismissing.(nt.seqprep) .& .!ismissing.(nt.eeg_age)) && return true
			# Keep any remaining subjects that have a seqprep *prior* to an EEG
			srt = sortperm([parse(Int, replace(v, "m"=>"")) for v in nt.visit])
			any(findfirst(!ismissing, nt.seqprep[srt]) .<= findfirst(!ismissing, nt.eeg_age[srt]))
		end))
	end
end

long_sub = let
	wide_sub = select(
		leftjoin(
			select(unstack(eegmbo, "subject", "visit", "eeg_age"),
				   "subject", "3m"=>"eeg_3m", "6m"=> "eeg_6m", "12m"=>"eeg_12m"),
			select(unstack(eegmbo, "subject", "visit", "age"),
				   "subject", "3m"=>"seqprep_3m", "6m"=> "seqprep_6m", "12m"=>"seqprep_12m"),
		on="subject"),
		"subject", r"3m", r"6m", r"12m"
	)

	long_sub = DataFrame()
	for row in eachrow(wide_sub), tp in tps
		stool_age = row["seqprep_$tp"]
		eeg_age = row["eeg_$tp"]
		push!(long_sub, (; subject=row.subject, timepoint=tp, stool_age, eeg_age); cols=:union)
	end

	@chain long_sub begin
		subset!(AsTable(["stool_age", "eeg_age"])=> ByRow(nt-> !all(ismissing, nt)))
		transform!(AsTable(["stool_age", "eeg_age"])=> ByRow(nt-> minimum(skipmissing(values(nt))))=> "minage") 
		sort!("minage")
	end
end


# Data for the EEG timeseries panel in Figure 1
# comes from separate analysis by Emma Margolis.

timeseries = mapreduce(vcat, ("3m", "6m", "12m")) do tp
    df = DataFrame(XLSX.readtable("data/alltimepoints_timeseries_figure.xlsx", "$tp Timeseries Calculation"; infer_eltypes=true))[:, 1:6]
    rename!(df, ["ms", "mean", "se", "lower", "upper", "std"])
    df.timepoint .= tp
    df[end, 1] = df[end-1, 1] + 1
    dropmissing!(df)
end

# Load gene functions files for samples that we have microbiomes for:


unirefs_by_sample = let ss = mapreduce(samplenames, vcat, (concurrent_3m, concurrent_6m, concurrent_12m))
    files = filter(f-> contains(basename(f), "genefamilies.tsv") && replace(basename(f), r"(SEQ\d+).+"=> s"\1") ∈ ss,
		readdir("/grace/sequencing/processed/mgx/humann/main/"; join=true)
    )
    mapreduce(vcat, files) do f
	df = CSV.read(f, DataFrame)
	rename!(df, ["feature", "abundance"])
	subset!(df, "feature"=> ByRow(f-> !contains(f, '|')))
	sample = replace(basename(f), r"(SEQ\d+).+"=> s"\1")
	df.sample .= sample
	df
    end
end

concurrent_unirefs = let 
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


set!(concurrent_unirefs, mapreduce(get, vcat, (concurrent_3m, concurrent_6m, concurrent_12m)))
unirefs_pco = pcoa(concurrent_unirefs)

concurrent_species = commjoin(concurrent_3m, concurrent_6m, concurrent_12m)
species_pco = pcoa(concurrent_species)


eeg = load_eeg()
eeg.peak_latency_P1_corrected = eeg.peak_latency_P1 .- eeg.peak_latency_N1
eeg.peak_latency_N2_corrected = eeg.peak_latency_N2 .- eeg.peak_latency_P1
eeg.peak_amp_P1_corrected = eeg.peak_amp_P1 .- eeg.peak_amp_N1
eeg.peak_amp_N2_corrected = eeg.peak_amp_N2 .- eeg.peak_amp_P1
rename!(eeg, "age"=> "eeg_age")
mbo = load_microbiome(eeg.subject)
transform!(mbo, "visit"=>ByRow(v-> replace(v, "mo"=>"m"))=> "visit")

eegmbo = @chain eeg begin
    select("subject", "timepoint"=>"visit", Cols(:))
    unique!(["subject", "timepoint"])
    outerjoin(mbo; on = ["subject", "visit"])
    groupby("subject")
    subset(AsTable(["visit", "seqprep", "eeg_age"])=> (nt-> begin
	# skip subjects without at least 1 eeg and at least 1 seqprep
	(all(ismissing, nt.seqprep) || all(ismissing(nt.eeg_age))) && return false
	# keep subjects with at least 1 concurrent eeg/seqprep
	any(.!ismissing.(nt.seqprep) .& .!ismissing.(nt.eeg_age)) && return true
	# Keep any remaining subjects that have a seqprep *prior* to an EEG
	srt = sortperm([parse(Int, replace(v, "m"=>"")) for v in nt.visit])
	any(findfirst(!ismissing, nt.seqprep[srt]) .<= findfirst(!ismissing, nt.eeg_age[srt]))
    end))
end

wide_sub = select(leftjoin(
    select(unstack(eegmbo, "subject", "visit", "eeg_age"), "subject", "3m"=>"eeg_3m", "6m"=> "eeg_6m", "12m"=>"eeg_12m"),
    select(unstack(eegmbo, "subject", "visit", "age"), "subject", "3m"=>"seqprep_3m", "6m"=> "seqprep_6m", "12m"=>"seqprep_12m"),
    on="subject"), "subject", r"3m", r"6m", r"12m"
)

# Load feature set enrichments

fsea_df = CSV.read("data/outputs/fsea/concurrent_consolidated_fsea.csv", DataFrame)
future6m_fsea_df = CSV.read("data/outputs/fsea/future6m_consolidated_fsea.csv", DataFrame)
future12m_fsea_df = CSV.read("data/outputs/fsea/future12m_consolidated_fsea.csv", DataFrame)

geneset_order = [
	"Acetate synthesis",
	"Acetate degradation",
	"Propionate synthesis",
	"Propionate degradation",
	"Butyrate synthesis",
	"Butyrate degradation",
	"GABA synthesis",
	"GABA degradation",
	"Glutamate synthesis",
	"Glutamate degradation",
	"Menaquinone synthesis",
	"Menaquinone degradation",
	"Inositol synthesis",
	"Inositol degradation",
	"Tryptophan synthesis",
	"Tryptophan degradation",
	"p-Cresol synthesis",
	"p-Cresol degradation",
	"Isovaleric acid synthesis",
	"Isovaleric acid degradation",
	"Quinolinic acid synthesis",
	"Quinolinic acid degradation",
	"S-Adenosylmethionine synthesis",
	"S-Adenosylmethionine degradation",
	"17-beta-Estradiol synthesis",
	"17-beta-Estradiol degradation",
	"DOPAC synthesis",
	"DOPAC degradation",
	"ClpB"
]
gssig = intersect(geneset_order,
	unique(subset(vcat(fsea_df, future6m_fsea_df, future12m_fsea_df), "q₀" => ByRow(<(0.2))).geneset)
)
gsidx = geneset_index(gssig,geneset_order)

alllms = groupby(mapreduce(vcat, Iterators.product(eeg_features, [tps..., ftps...])) do (feat, tp)
	df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
	rename!(df, "Name"=>"eeg_feature")
	df.timepoint .= tp
	df
end, ["eeg_feature", "timepoint"])

for gs in gssig
	@info gs
	transform!(alllms, "feature"=> ByRow(u-> replace(u, "UniRef90_"=>"") ∈ na_map[gs])=> gs; ungroup=false)
end
##

topspecies = CSV.read("data/outputs/fsea_top_species.csv", DataFrame)

humann_files = let
	seqs = Set(mapreduce(samplenames, vcat,
		(concurrent_3m, concurrent_6m, concurrent_12m,
		 future_3m6m, future_3m12m, future_6m12m)))
	na_unirefs = Set(reduce(union, values(na_map)))
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
	
	for gs in gssig
		@info gs
		transform!(humann_files, "uniref"=> ByRow(u-> u ∈ na_map[gs])=> gs)
	end
	groupby(humann_files, "species")
end

topspec_idx = Dict(s=> i for (i, s) in enumerate(sort(unique(filter(!=("other"), topspecies.species)))))
topspec_unirefs_idx = Dict(u=> i for (i, u) in enumerate(unique(select(humann_files, "uniref").uniref)))

topspec_eeg_cormat_3m = zeros(length(keys(topspec_idx)), 6)
topspec_eeg_cormat_6m = zeros(length(keys(topspec_idx)), 6)
topspec_eeg_cormat_12m = zeros(length(keys(topspec_idx)), 6)

for (i, feat) in enumerate(eeg_features), spc in keys(topspec_idx)
	if haskey(concurrent_3m.fidx, spc)
		topspec_eeg_cormat_3m[topspec_idx[spc], i] = cor(
			vec(abundances(concurrent_3m[spc, :])),
			get(concurrent_3m, Symbol(feat))
		)
	end
	if haskey(concurrent_6m.fidx, spc)
		topspec_eeg_cormat_6m[topspec_idx[spc], i] = cor(
			vec(abundances(concurrent_6m[spc, :])),
			get(concurrent_6m, Symbol(feat))
		)
	end
	if haskey(concurrent_12m.fidx, spc)
		topspec_eeg_cormat_12m[topspec_idx[spc], i] = cor(
			vec(abundances(concurrent_12m[spc, :])),
			get(concurrent_12m, Symbol(feat))
		)
	end
end

topspec_eeg_cormat_full = hcat(topspec_eeg_cormat_3m, topspec_eeg_cormat_6m, topspec_eeg_cormat_12m)

topspec_eeg_hcl_3m_row = hclust(
	pairwise(Euclidean(), topspec_eeg_cormat_3m; dims=1);
	linkage=:complete,
	branchorder=:optimal
)
topspec_eeg_hcl_3m_col = hclust(
	pairwise(Euclidean(), topspec_eeg_cormat_3m; dims=2);
	linkage=:complete,
	branchorder=:optimal
)

topspec_eeg_hcl_6m_row = hclust(
	pairwise(Euclidean(), topspec_eeg_cormat_6m; dims=1);
	linkage=:complete,
	branchorder=:optimal
)
topspec_eeg_hcl_6m_col = hclust(
	pairwise(Euclidean(), topspec_eeg_cormat_6m; dims=2);
	linkage=:complete,
	branchorder=:optimal
)

topspec_eeg_hcl_12m_row = hclust(
	pairwise(Euclidean(), topspec_eeg_cormat_12m; dims=1);
	linkage=:complete,
	branchorder=:optimal
)
topspec_eeg_hcl_12m_col = hclust(
	pairwise(Euclidean(), topspec_eeg_cormat_12m; dims=2);
	linkage=:complete,
	branchorder=:optimal
)

topspec_eeg_hcl_full_row = hclust(
    pairwise(Euclidean(), topspec_eeg_cormat_full; dims=2);
	linkage=:complete,
	branchorder=:optimal
)
topspec_eeg_hcl_full_col = hclust(
    pairwise(Euclidean(), topspec_eeg_cormat_full; dims=1);
	linkage=:complete,
	branchorder=:optimal
)

##

topspec_unirefs_abmat = spzeros(length(keys(topspec_unirefs_idx)), length(topspec_idx))

for subdf in humann_files
	spec = first(subdf.species)
	haskey(topspec_idx, spec) || continue
	abs = combine(groupby(subdf, "uniref"), "uniref" => first => "uniref", "abundance"=> mean => "abmean")
	for row in eachrow(abs)
		topspec_unirefs_abmat[topspec_unirefs_idx[row.uniref], topspec_idx[spec]] = row.abmean
	end
end
 
topspec_unirefs_hcl_row = hclust(
    pairwise(Euclidean(), topspec_unirefs_abmat; dims=2);
	linkage=:complete,
	branchorder=:optimal
)
topspec_unirefs_hcl_col = hclust(
    pairwise(Euclidean(), topspec_unirefs_abmat; dims=1);
	linkage=:complete,
	branchorder=:optimal
)

##



# ## Plotting
# 
# Each figure plot has a common structure.
#
# 1. Initialize the figure with `Figure()`
# 2. Lay out the figure using `GridLayout()` to get the overall structure.
#    Here, I'm using variable names describing content rather than panels
#    in case things need to change later.
# 3. Create and plot into specific axes (`Axis()`), add legends etc.
# 4. Save figures to SVG and PNG formats.
#
# ### Colors

colormap_age = :viridis
colors_timepoints = [tp=> c for (tp, c) in zip(tps, cgrad(colormap_age)[[0., 0.4, 0.8]])]
colors_sampletype = [st=> c for (st, c) in zip(["stool", "eeg", "both"], cgrad(:tab10; categorical=true)[[1,4,5]])]
colors_eeg_kind = [kd=> c for (kd, c) in zip(["amp", "latency"], cgrad(:tab10; categorical=true)[[2,3]])]
colors_eeg_peaks = [pk=> c for (pk, c) in zip(["N1", "P1c", "N2c"], cgrad(:Accent_4; categorical=true)[2:4])]

colormap_sig = :PuOr
colors_sig = cgrad(colormap_sig, 11; rev = true, categorical=true)[[1,3,4,8,9,11]]
insert!(colors_sig, 4, colorant"white")




# ### Main figures
#
# #### Figure 1

figure1 = Figure(; size = (1200,800))

grid_cohort = GridLayout(figure1[1,1:2])
grid_longsamples = GridLayout(figure1[2:3, 1])
grid_eeg_curves = GridLayout(figure1[3,2])

grid_pcoas = GridLayout(figure1[2, 2])

ax_cohort = Axis(grid_cohort[1,1]; aspect = DataAspect(), alignmode=Outside())
# Legend(grid_cohort[2,1],
# 	[MarkerElement(; marker=:rect, color=c[2]) for c in colors_timepoints],
# 	["visit 1", "visit 2", "visit 3"];
# 	orientation=:horizontal, tellwidth=false, tellheight=true
# )

ax_eeg_curves = Axis(grid_eeg_curves[1,1]; xlabel="time (ms) relative to stimulus onset",
		     ylabel="voltage (μV)"
)
ax_pcoa_spec = Axis(grid_pcoas[1,1])
ax_pcoa_func = Axis(grid_pcoas[1,2])

ax_longsamples = Axis(grid_longsamples[3,1]; xlabel="age (months)", ylabel = "subject")
ax_eeg_hist = Axis(grid_longsamples[1,1]; ylabel="density", title="eeg")
ax_stool_hist = Axis(grid_longsamples[2,1]; ylabel="density", title="stool")


# ##### Cohort cartoon
#
# This image was made externally.
# Note: to get this into the right position/size
# without distortion, edits have to be made to the original image.

image!(ax_cohort, rotr90(load("manuscript/mainfigures/figure_draft1.png")))

hidedecorations!(ax_cohort)
hidespines!(ax_cohort)

# ##### EEG longitudinal curves


datmean = data(timeseries) * mapping(:ms, :mean, color=:timepoint)
datlower = data(timeseries) * mapping(:ms, :lower, color=:timepoint)
datupper = data(timeseries) * mapping(:ms, :upper, color=:timepoint)

let datage = data(subset(eegmbo, "eeg_age"=>ByRow(!ismissing))) * mapping(:eeg_age=> "age (months)", color="visit")
	draw!(ax_eeg_hist, datage * AlgebraOfGraphics.density(); palettes=(; color=colors_timepoints))
end

let datage = data(subset(eegmbo, "age"=>ByRow(!ismissing))) * mapping(:age=> "age (months)", color="visit")
	draw!(ax_stool_hist, datage * AlgebraOfGraphics.density(); palettes=(; color=colors_timepoints))
end

mpl = draw!(ax_eeg_curves, datmean * visual(Lines); palettes=(; color=colors_timepoints))
dpl = draw!(ax_eeg_curves, datlower * visual(Lines; linestyle=:dash); palettes=(; color=colors_timepoints))
draw!(ax_eeg_curves, datupper * visual(Lines; linestyle=:dash); palettes=(; color=colors_timepoints))

Legend(grid_eeg_curves[1,1],
     [LineElement(; color=:gray), LineElement(; color=:gray, linestyle=:dash)],
     ["mean", "+/- S.E."];
	 tellwidth=false, halign=:right, valign=:top, margin=(10,10,10,10)
) 

# ##### PCoAs


plt_pcoa_spec = plot_pcoa!(ax_pcoa_spec, species_pco; color=get(concurrent_species, :age), colormap=colormap_age)
plt_pcoa_func = plot_pcoa!(ax_pcoa_func, unirefs_pco; color=get(concurrent_unirefs, :age), colormap=colormap_age)

Colorbar(grid_pcoas[1,3];
	limits=extrema(skipmissing(eegmbo.age)), label = "Age (months)", colormap=colormap_age,
)

# ##### Longitudinal


subind = Dict(s=> i for (i,s) in enumerate(unique(long_sub.subject)))

samples_colors = [s => c for (s, c) in zip(
	("stool", "eeg", "both"),
	(colorant"dodgerblue", colorant"firebrick", colorant"slateblue"))
]

let gdf = groupby(long_sub, "subject")
	for k in keys(gdf)
		stools = gdf[k].stool_age
		eegs = gdf[k].eeg_age

		y = subind[k.subject]
		xs = map(zip(stools, eegs)) do (s,e)
			ismissing(s) && return e
			ismissing(e) && return s
			return mean([e,s])
		end
		cs = map(zip(stools, eegs)) do (s,e)
			ismissing(s) && return colors_sampletype[2][2]
			ismissing(e) && return colors_sampletype[1][2]
			return colors_sampletype[3][2]
		end
		scatter!(ax_longsamples, xs, fill(y, length(xs)); color=cs)
	 
	 
		lines!(ax_longsamples, [extrema(skipmissing([stools;eegs]))...], [y,y]; linestyle=:dash, color = :gray)
		for s in filter(!ismissing, stools)
			lines!(ax_longsamples, [s,s], [y+0.4, y-0.4]; color=colors_sampletype[1][2])
		end
		for e in filter(!ismissing, eegs)
			lines!(ax_longsamples, [e,e], [y+0.4, y-0.4]; color=colors_sampletype[2][2])
		end
	end
end

# Legend(grid_longsamples[3,1],
# 	[MarkerElement(; marker=:circle, color=colors_sampletype[i][2]) for i in 1:3],
# 	["stool", "eeg", "both"];
# 	tellwidth=false, valign=:top, halign=:left, margin=(10,10,10,10)
# )

ylims!(ax_longsamples, -2, length(unique(long_sub.subject)) + 2)


# tweak some visuals

rowsize!(grid_longsamples, 3, Relative(5/6))
colsize!(figure1.layout, 1, Relative(1/3))

hidedecorations!(ax_eeg_curves, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_pcoa_spec, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_pcoa_func, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_longsamples, ticks=false, ticklabels=false, label=false)
hidexdecorations!.((ax_stool_hist, ax_eeg_hist))
hideydecorations!.((ax_stool_hist, ax_eeg_hist); ticks=false, ticklabels=false, label=false)
linkyaxes!(ax_eeg_hist, ax_stool_hist)

save("/home/kevin/Downloads/figure1-inprogress.png", figure1)

# ##### Figure 2


figure2 = Figure(; size = (1100,600))

grid_fsea_dots = GridLayout(figure2[1,1])
# grid_fsea_heatmaps = GridLayout(figure2[1,2])
# grid_bug_heatmaps = GridLayout(figure2[2, 1])

gs_interval = 6
tp_interval = 1.5
dotsxlim = 3
fsea_marker_size=6

grid_fsea_latency=GridLayout(grid_fsea_dots[1,1])
ax_lat = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_fsea_latency[1,i];
				 xlabel = "z",
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 )
		r = Axis(grid_fsea_latency[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(string.(3:-1:1); outer=length(tickrange)))
				 )
		(l, r)
	end
end

grid_fsea_amplitude=GridLayout(grid_fsea_dots[1,2])
ax_amp = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_fsea_amplitude[1,i];
				 xlabel = "z",
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 )
		r = Axis(grid_fsea_amplitude[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(string.(3:-1:1); outer=length(tickrange)))
				 )
		(l, r)
	end
end

for axs in (ax_lat, ax_amp)
	for (i, (axl, axr)) in enumerate(axs)
		ylims!.((axl,axr), gs_interval - 2*tp_interval, gs_interval * length(gssig) + 2*tp_interval)
		xlims!(axl, -dotsxlim, dotsxlim)
		hidexdecorations!(axl; ticks=false, ticklabels=false, label=false)
		hideydecorations!(axl, ticks=i != 1, ticklabels=i != 1)
		hideydecorations!(axr, ticks=i != 3, ticklabels=i != 3)
		hidexdecorations!(axr)

		for gs in gssig
			mid = gsidx[gs] * gs_interval
			isodd(gsidx[gs]) || continue
			span = gs_interval / 2
			poly!(axl, Point2f.([(-dotsxlim, mid - span),
								 (dotsxlim, mid - span),
								 (dotsxlim, mid + span),
								 (-dotsxlim, mid + span)]);
				  color=("gray80", 0.4))
		end

	end
end
hideydecorations!(ax_amp[1][1], ticks=false)

Legend(grid_fsea_dots[2,1:2], 
	   [[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 1:3],
		[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 5:7]],
	   [["q < 0.01", "q < 0.1", "q < 0.2"], 
		["q < 0.2", "q < 0.1", "q < 0.01"]],
	   ["(-)", "(+)"];
	   orientation=:horizontal, tellheight=true, tellwidth=false)
Label(grid_fsea_dots[0,1], "Latency"; fontsize=20, tellwidth=false)
Label(grid_fsea_dots[0,2], "Amplitude"; fontsize=20, tellwidth=false)

#-


for (j, feat) in enumerate(filter(contains("latency"), eeg_features))
	gdf = groupby(fsea_df, "eeg_feature")
    featstr = replace(feat, "peak_latency_"=>"", "_corrected"=>"c")
    @warn "$featstr"
    for (i, tp) in enumerate(tps)
        @info "$tp"

		df =  alllms[(; eeg_feature=feat, timepoint=tp)]
		allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        for row in eachrow(subdf)
			yidx = df[!, row.geneset]
			xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
			ys = df.z[yidx]
			xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_lat[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_lat[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_lat[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
end


for (j, feat) in enumerate(filter(contains("amp"), eeg_features))
	gdf = groupby(fsea_df, "eeg_feature")
    featstr = replace(feat, "peak_amp"=>"", "_corrected"=>"c")
    @warn "$featstr"
    for (i, tp) in enumerate(tps)
        @info "$tp"

		df =  alllms[(; eeg_feature=feat, timepoint=tp)]
		allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        # lines!(ax_lat, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        for row in eachrow(subdf)
			# haskey(gsidx, row.geneset) || continue
			yidx = df[!, row.geneset]
			xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
			ys = df.z[yidx]
			xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_amp[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_amp[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_amp[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
end

#-


# ax_bug_heatmap = Axis(grid_bug_heatmaps[1,1]; 
# 	yticks = (1:size(topspec_eeg_cormat_full, 1),
# 			  collect(keys(topspec_idx))[sortperm(collect(values(topspec_idx)))][topspec_eeg_hcl_full_col.order]
# 	),
# )
# hidexdecorations!(ax_bug_heatmap)
# ax_bug_annotations = Axis(grid_bug_heatmaps[2,1])
# hidedecorations!(ax_bug_annotations)
#
# Legend(grid_bug_heatmaps[1,2],
# 	   [MarkerElement(; color=c, marker=:rect) for c in mapreduce(
#			x-> x[2], vcat, [colors_timepoints; colors_eeg_peaks; colors_eeg_kind]
#			)
#	   ],
# 	   mapreduce(x-> x[1], vcat, [colors_timepoints; colors_eeg_peaks; colors_eeg_kind])
# )
#
# # ax_bug_heatmap_6m = Axis(grid_bug_heatmaps[1,2]; title="6m", )
# # ax_bug_heatmap_12m = Axis(grid_bug_heatmaps[1,3]; title="12m")
#
# heatmap!(ax_bug_heatmap,
#		collect(topspec_unirefs_abmat)[topspec_unirefs_hcl_col.order, topspec_unirefs_hcl_row.order]
# )
#
# for (i, (feat, tp)) in enumerate(
# 	vec((collect(Iterators.product(eeg_features, tps))))[topspec_eeg_hcl_full_row.order]
# )
# 	_, kind, peak = split(replace(feat, "_corrected"=>"c"), "_")
# 	kindcol = Dict(colors_eeg_kind)[kind]
# 	peakcol = Dict(colors_eeg_peaks)[peak]
# 	tpcol = Dict(colors_timepoints)[tp]
# 	poly!(ax_bug_annotations, Point2f.([
# 				(i-0.5, 0), (i+0.5, 0), (i+0.5, 1), (i-0.5, 1)
# 			]); color=kindcol)
#
# 	poly!(ax_bug_annotations, Point2f.([
# 				(i-0.5, 1), (i+0.5, 1), (i+0.5, 2), (i-0.5, 2)
# 			]); color=peakcol)
#
# 	poly!(ax_bug_annotations, Point2f.([
# 				(i-0.5, 2), (i+0.5, 2), (i+0.5, 3), (i-0.5, 3)
# 			]); color=tpcol)
#
# end
# tightlimits!(ax_bug_annotations)
# rowsize!(figure2.layout, 2, Relative(1/3))
# colgap!.((grid_fsea_latency, grid_fsea_amplitude), 2)
#
# for feat in filter(contains("latency"), feats)
save("/home/kevin/Downloads/figure2-inprogress.png", figure2)

############
# Figure 3 #
############

figure3 = Figure(; size = (1050, 750))


grid_future_violins = GridLayout(figure3[1,1])
grid_fsea_dots = GridLayout(figure3[2,1])

gs_interval = 6
tp_interval = 1.5
dotsxlim = 3
fsea_marker_size=6

grid_fsea_latency=GridLayout(grid_fsea_dots[1,1])
ax_lat = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_fsea_latency[1,i];
				 xlabel = "z",
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 )
		r = Axis(grid_fsea_latency[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(string.(3:-1:1); outer=length(tickrange)))
				 )
		(l, r)
	end
end

grid_fsea_amplitude=GridLayout(grid_fsea_dots[1,2])
ax_amp = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_fsea_amplitude[1,i];
				 xlabel = "z",
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 )
		r = Axis(grid_fsea_amplitude[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(string.(3:-1:1); outer=length(tickrange)))
				 )
		(l, r)
	end
end

for axs in (ax_lat, ax_amp)
	for (i, (axl, axr)) in enumerate(axs)
		ylims!.((axl,axr), gs_interval - 2*tp_interval, gs_interval * length(gssig) + 2*tp_interval)
		xlims!(axl, -dotsxlim, dotsxlim)
		hidexdecorations!(axl; ticks=false, ticklabels=false, label=false)
		hideydecorations!(axl, ticks=i != 1, ticklabels=i != 1)
		hideydecorations!(axr, ticks=i != 3, ticklabels=i != 3)
		hidexdecorations!(axr)

		for gs in gssig
			mid = gsidx[gs] * gs_interval
			isodd(gsidx[gs]) || continue
			span = gs_interval / 2
			poly!(axl, Point2f.([(-dotsxlim, mid - span),
								 (dotsxlim, mid - span),
								 (dotsxlim, mid + span),
								 (-dotsxlim, mid + span)]);
				  color=("gray80", 0.4))
		end

	end
end
hideydecorations!(ax_amp[1][1], ticks=false)

Legend(grid_fsea_dots[2,1:2], 
	   [[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 1:3],
		[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 5:7]],
	   [["q < 0.01", "q < 0.1", "q < 0.2"], 
		["q < 0.2", "q < 0.1", "q < 0.01"]],
	   ["(-)", "(+)"];
	   orientation=:horizontal, tellheight=true, tellwidth=false)
Label(grid_fsea_dots[0,1], "Latency"; fontsize=20, tellwidth=false)
Label(grid_fsea_dots[0,2], "Amplitude"; fontsize=20, tellwidth=false)

#-


for (j, feat) in enumerate(filter(contains("latency"), eeg_features))
	gdf = groupby(fsea_df, "eeg_feature")
    featstr = replace(feat, "peak_latency_"=>"", "_corrected"=>"c")
    @warn "$featstr"
    for (i, tp) in enumerate(tps)
        @info "$tp"

		df =  alllms[(; eeg_feature=feat, timepoint=tp)]
		allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        for row in eachrow(subdf)
			yidx = df[!, row.geneset]
			xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
			ys = df.z[yidx]
			xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_lat[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_lat[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_lat[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
end


for (j, feat) in enumerate(filter(contains("amp"), eeg_features))
	gdf = groupby(fsea_df, "eeg_feature")
    featstr = replace(feat, "peak_amp"=>"", "_corrected"=>"c")
    @warn "$featstr"
    for (i, tp) in enumerate(tps)
        @info "$tp"

		df =  alllms[(; eeg_feature=feat, timepoint=tp)]
		allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        # lines!(ax_lat, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        for row in eachrow(subdf)
			# haskey(gsidx, row.geneset) || continue
			yidx = df[!, row.geneset]
			xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
			ys = df.z[yidx]
			xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_amp[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_amp[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_amp[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
end

#-


