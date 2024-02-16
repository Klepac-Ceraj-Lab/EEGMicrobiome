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
using ThreadsX
using Chain
using XLSX
using FileIO
using CSV
using DataFrames
using Distributions
using Microbiome
using BiobakeryUtils
using CairoMakie
using AlgebraOfGraphics
using SparseArrays

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
# ### Main figures
#
# #### Figure 1

figure1 = Figure(; size = (1200,800))

grid_cohort = GridLayout(figure1[1,1])
grid_eeg_curves = GridLayout(figure1[2,1])

grid_pcoas = GridLayout(figure1[1:2, 2])
grid_longsamples = GridLayout(figure1[1:2, 3])

ax_cohort = Axis(grid_cohort[1,1]; aspect = DataAspect(), alignmode=Outside())

ax_eeg_curves = Axis(grid_eeg_curves[1,1]; xlabel="time (ms) relative to stimulus onset",
		     ylabel="voltage (μV)"
)
ax_pcoa_spec = Axis(grid_pcoas[1,1])
ax_pcoa_func = Axis(grid_pcoas[2,1])

ax_longsamples = Axis(grid_longsamples[3:5,1]; xlabel="age (months)", ylabel = "subject")
ax_vep_hist = Axis(grid_longsamples[1,1]; ylabel="samples (fraction)")
ax_stool_hist = Axis(grid_longsamples[2,1]; ylabel="samples (fraction)")


# ##### Cohort cartoon
#
# This image was made externally.
# Note: to get this into the right position/size
# without distortion, edits have to be made to the original image.

image!(ax_cohort, rotr90(load("manuscript/mainfigures/eegfig1a.png")))

hidedecorations!(ax_cohort)
hidespines!(ax_cohort)

# ##### EEG longitudinal curves


datmean = data(timeseries) * mapping(:ms, :mean, color=:timepoint)
datlower = data(timeseries) * mapping(:ms, :lower, color=:timepoint)
datupper = data(timeseries) * mapping(:ms, :upper, color=:timepoint)

eeg_timepoint_colors = [tp=> c for (tp, c) in zip(("3m", "6m", "12m"), cgrad(:solar, 3; categorical=:true))]
stoo_timepoint_colors = [tp=> c for (tp, c) in zip(("3m", "6m", "12m"), cgrad(:deep, 3; categorical=:true))]

let datage = data(subset(eegmbo, "eeg_age"=>ByRow(!ismissing))) * mapping(:eeg_age=> "age (months)", color="visit")
	draw!(ax_eeg_hist, datage * AlgebraOfGraphics.density(); palettes=(; color=eeg_timepoint_colors))
end

let datage = data(subset(eegmbo, "age"=>ByRow(!ismissing))) * mapping(:age=> "age (months)", color="visit")
	draw!(ax_stool_hist, datage * AlgebraOfGraphics.density(); palettes=(; color=stool_timepoint_colors))
end

mpl = draw!(ax_eeg_curves, datmean * visual(Lines); palettes=(; color=eeg_timepoint_colors))
dpl = draw!(ax_eeg_curves, datlower * visual(Lines; linestyle=:dash); palettes=(; color=eeg_timepoint_colors))
draw!(ax_eeg_curves, datupper * visual(Lines; linestyle=:dash); palettes=(; color=eeg_timepoint_colors))

Legend(grid_eeg_curves[2,1],
    [[LineElement(; color=timepoint_colors[i][2]) for i in 1:3],
     [LineElement(; color=:gray), LineElement(; color=:gray, linestyle=:dash)]
    ],
    [["visit 1", "visit 2", "visit 3"],
     ["mean", "+/- S.E."]],
    ["visit", "value"];
    orientation=:horizontal, tellheight=true, tellwidth=false,
) 

# ##### PCoAs


plt_pcoa_spec = plot_pcoa!(ax_pcoa_spec, species_pco; color=get(concurrent_species, :age))
Colorbar(grid_pcoas[1,2], plt_pcoa_spec; label = "age (months)")


plt_pcoa_func = plot_pcoa!(ax_pcoa_func, unirefs_pco; color=get(concurrent_unirefs, :age))
Colorbar(grid_pcoas[2,2], plt_pcoa_func; label="age (months)")

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
			ismissing(s) && return samples_colors[2][2]
			ismissing(e) && return samples_colors[1][2]
			return samples_colors[3][2]
		end
		scatter!(ax_longsamples, xs, fill(y, length(xs)); color=cs)
	 
	 
		lines!(ax_longsamples, [extrema(skipmissing([stools;eegs]))...], [y,y]; linestyle=:dash, color = :gray)
		for s in filter(!ismissing, stools)
			lines!(ax_longsamples, [s,s], [y+0.4, y-0.4]; color=samples_colors[1][2])
		end
		for e in filter(!ismissing, eegs)
			lines!(ax_longsamples, [e,e], [y+0.4, y-0.4]; color=samples_colors[2][2])
		end
	end
end

Legend(grid_longsamples[1,1],
	[MarkerElement(; marker=:circle, color=c[2]) for c in samples_colors],
	[c[1] for c in samples_colors];
	halign=:left, valign=:top, tellwidth=false, margin=(10,10,10,10)
)


ylims!(ax_longsamples, -2, length(unique(long_sub.subject)) + 2)


# tweak some visuals

hidedecorations!(ax_eeg_curves, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_pcoa_spec, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_pcoa_func, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_longsamples, ticks=false, ticklabels=false, label=false)

save("/home/kevin/Downloads/figure1-inprogress.png", figure1)
