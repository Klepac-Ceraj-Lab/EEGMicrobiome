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
using DataFrames
using Distributions
using Microbiome
using BiobakeryUtils
using CairoMakie
using AlgebraOfGraphics
using AlgebraOfGraphics: categorical
using SparseArrays
using Clustering
using Distances

# Next, we'll load the data.

tps = ("v1","v2", "v3")
ftps = ("v1v2", "v1v3", "v2v3")

mdata = load_cohorts()

long_sub = let
	wide_sub = select(
		leftjoin(
			select(unstack(mdata, "subject_id", "visit", "eeg_age"),
				   "subject_id", "v1"=>"eeg_v1", "v2"=> "eeg_v2", "v3"=>"eeg_v3"),
			select(unstack(mdata, "subject_id", "visit", "stool_age"),
				   "subject_id", "v1"=>"seqprep_v1", "v2"=> "seqprep_v2", "v3"=>"seqprep_v3"),
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
		subset!(AsTable(["stool_age", "eeg_age"])=> ByRow(nt-> !all(ismissing, nt)))
		transform!(AsTable(["stool_age", "eeg_age"])=> ByRow(nt-> minimum(skipmissing(values(nt))))=> "minage") 
		sort!("minage")
	end
end


v1 = get_cohort(mdata, "v1")
v2 = get_cohort(mdata, "v2")
v3 = get_cohort(mdata, "v3")
v1v2 = get_cohort(mdata, "v1v2")
v1v3 = get_cohort(mdata, "v1v3")
v2v3 = get_cohort(mdata, "v2v3")

##

eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency"=>"amp")]

na_map = FeatureSetEnrichments.get_neuroactive_unirefs()

# Data for the EEG timeseries panel in Figure 1
# comes from separate analysis by Emma Margolis.

vep_timeseries = mapreduce(vcat, enumerate(("3m", "6m", "12m"))) do (i, tp)
    df = DataFrame(XLSX.readtable("data/alltimepoints_timeseries_figure.xlsx", "$tp Timeseries Calculation"; infer_eltypes=true))[:, 1:6]
    tpv = "v$i"
    rename!(df, ["ms", "mean", "se", "lower", "upper", "std"])
    df.timepoint .= tpv
    df[end, 1] = df[end-1, 1] + 1
    dropmissing!(df)
end

# Load gene functions files for samples that we have microbiomes for:

taxprofiles = let mpa = metaphlan_profiles(String.(skipmissing(mdata.taxprofile)), :species)
	CommunityProfile(abundances(mpa), features(mpa), MicrobiomeSample.(replace.(samplenames(mpa), r"_S\d++_profile"=> "")))
end

set!(taxprofiles, select(mdata, "seqprep"=> "sample", Cols(:)))

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

unirefs_pco = pcoa(unirefs)
species_pco = pcoa(taxprofiles)

fsea_df = CSV.read("data/outputs/fsea/concurrent_consolidated_fsea.csv", DataFrame)
futfsea_df = CSV.read("data/outputs/fsea/future_consolidated_fsea.csv", DataFrame)

geneset_types = (;
	neurotransmitters = [
		"GABA synthesis",
		"GABA degradation",
		"Glutamate synthesis",
		"Glutamate degradation",
		],
	nt_metabolism = [
		# tryptophan/serotonin stuff
		"Tryptophan synthesis",
		"Tryptophan degradation",
		"Quinolinic acid synthesis",
		"Quinolinic acid degradation",
		# dopamine stuff
		"DOPAC synthesis",
		"DOPAC degradation",
		],
	scfa = [
		"Acetate synthesis",
		"Acetate degradation",
		"Propionate synthesis",
		"Propionate degradation",
		"Butyrate synthesis",
		"Butyrate degradation",
		"Isovaleric acid synthesis",
		"Isovaleric acid degradation",
		],
	other = [
		"Menaquinone synthesis",
		"Menaquinone degradation",
		"Inositol synthesis",
		"Inositol degradation",
		"p-Cresol synthesis",
		"p-Cresol degradation",
		"S-Adenosylmethionine synthesis",
		"S-Adenosylmethionine degradation",
		"17-beta-Estradiol synthesis",
		"17-beta-Estradiol degradation",
		"ClpB"
		]
)
gs_types_rev = Dict(gs=>t for t in keys(geneset_types) for gs in geneset_types[t]) 

geneset_order = vcat(geneset_types...)
gssig = intersect(geneset_order,
	unique(subset(vcat(fsea_df, futfsea_df), "q₀" => ByRow(<(0.2))).geneset)
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

colors_gstypes = Dict(gst => c for (gst, c) in zip(keys(geneset_types), cgrad(:tab10; categorical=true)[[10,3,2,9]]))

# ### Tables
#
# #### Table 2

table2 = subset(fsea_df, "timepoint"=> ByRow(t-> t ∈ tps),
		   "q₀"=> ByRow(<(0.2)),
		   "geneset"=> ByRow(g-> g ∈ geneset_order)
)
transform!(table2, "eeg_feature" => ByRow(f-> begin
			(_, p, n...) = split(f, "_")
			n = first(n)
			return (; feature=p, peak=n)
	end) => ["feature", "peak"]
)
transform!(table2, "timepoint"=> ByRow(x-> Dict("v1"=> 1, "v2"=> 2, "v3"=> 3)[x]) => "visit")
select!(table2, "visit", "geneset", "feature", "peak", "es"=> "E.S.", "q₀"=> "q value")
sort!(table2, [
	"visit",
	order("geneset"; by = x-> findfirst(==(x), geneset_order)),
	"feature",
	# order("peak"; by = x-> findfirst(==(x), ("N1", "P1", "N2")))
	"E.S."
	]
)

for v in 1:3
	CSV.write("/home/kevin/Downloads/Table2_v$v.csv", subset(table2, "visit"=> ByRow(==(v))))
end


# #### Table 3

table3 = subset(futfsea_df,
		    "timepoint"=> ByRow(t-> t ∈ ftps),
		    "q₀"=> ByRow(<(0.2)),
		    "geneset"=> ByRow(g-> g ∈ geneset_order)
)
transform!(table3, "eeg_feature" => ByRow(f-> begin
			(_, p, n...) = split(f, "_")
			n = first(n)
			return (; feature=p, peak=n)
	end) => ["feature", "peak"]
)
transform!(table3, "timepoint"=> ByRow(x-> begin
		(st, vep) = match(r"^v(\d+)v(\d+)$", x).captures
		return (; stool = parse(Int, st), VEP = parse(Int, vep))
	end) => ["stool", "VEP"]
)
select!(table3, "stool", "VEP", "geneset", "feature", "peak", "es"=> "E.S.", "q₀"=> "q value")
sort!(table3, [
	"stool",
	"VEP",
	order("geneset"; by = x-> findfirst(==(x), geneset_order)),
	"feature",
	# order("peak"; by = x-> findfirst(==(x), ("N1", "P1", "N2")))
	"E.S."
	]
)
for (s, v) in zip((1,1,2), (2,3,3))
	CSV.write("/home/kevin/Downloads/Table3_$(s)_$(v).csv", subset(table3, "stool"=>ByRow(==(s)), "VEP"=>ByRow(==(v))))
end

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
		     ylabel="voltage (μV)",
			 yminorticksvisible=true,
			 yminorticks=IntervalsBetween(5)
)
ax_pcoa_spec = Axis(grid_pcoas[1,1])
ax_pcoa_func = Axis(grid_pcoas[1,2])

ax_longsamples = Axis(grid_longsamples[3,1];
					 xlabel ="age (months)",
					 xticks = 0:3:18,
					 limits = ((0., 19.), nothing),
					 ylabel = "subject"
				)
ax_eeg_hist = Axis(grid_longsamples[1,1]; ylabel="density", title="eeg")
ax_stool_hist = Axis(grid_longsamples[2,1]; ylabel="density", title="stool")


# ##### Cohort cartoon
#
# This image was made externally.
# Note: to get this into the right position/size
# without distortion, edits have to be made to the original image.

image!(ax_cohort, rotr90(load("manuscript/mainfigures/FinalConceptFigure.png")))

hidedecorations!(ax_cohort)
hidespines!(ax_cohort)

# ##### EEG longitudinal curves


datmean = data(vep_timeseries) * mapping(:ms, :mean, color=:timepoint)
datlower = data(vep_timeseries) * mapping(:ms, :lower, color=:timepoint)
datupper = data(vep_timeseries) * mapping(:ms, :upper, color=:timepoint)

let datage = data(subset(mdata, "eeg_age"=>ByRow(!ismissing))) * mapping(:eeg_age=> "age (months)", color="visit")
	draw!(ax_eeg_hist, datage * AlgebraOfGraphics.density(); palettes=(; color=colors_timepoints))
end

let datage = data(subset(mdata, "stool_age"=>ByRow(!ismissing))) * mapping(:stool_age=> "age (months)", color="visit")
	draw!(ax_stool_hist, datage * AlgebraOfGraphics.density(); palettes=(; color=colors_timepoints))
end

mpl = draw!(ax_eeg_curves, datmean * visual(Lines); palettes=(; color=colors_timepoints))
dpl = draw!(ax_eeg_curves, datlower * visual(Lines; linestyle=:dash); palettes=(; color=colors_timepoints))
draw!(ax_eeg_curves, datupper * visual(Lines; linestyle=:dash); palettes=(; color=colors_timepoints))

Legend(grid_eeg_curves[1,1],
     [LineElement(; color=:gray), LineElement(; color=:gray, linestyle=:dash)],
     ["mean", "+/- S.E."];
	 tellwidth=false, halign=:right, valign=:bottom, margin=(10,10,10,10)
) 

Legend(grid_eeg_curves[1,1],
    [LineElement(; color=c[2]) for c in colors_timepoints],
    ["Visit $i" for i in 1:3];
	tellwidth=false, halign=:right, valign=:top, margin=(10,10,10,10)
) 
# ##### PCoAs


plt_pcoa_spec = plot_pcoa!(ax_pcoa_spec, species_pco; color=get(taxprofiles, :stool_age), colormap=colormap_age)
plt_pcoa_func = plot_pcoa!(ax_pcoa_func, unirefs_pco; color=get(unirefs, :stool_age), colormap=colormap_age)

Colorbar(grid_pcoas[1,3];
	limits=extrema(skipmissing(mdata.stool_age)), label = "Age (months)", colormap=colormap_age,
)

# ##### Longitudinal


subind = Dict(s=> i for (i,s) in enumerate(unique(long_sub.subject_id)))

samples_colors = [s => c for (s, c) in zip(
	("stool", "eeg", "both"),
	(colorant"dodgerblue", colorant"firebrick", colorant"slateblue"))
]

let gdf = groupby(long_sub, "subject_id")
	for k in keys(gdf)
		stools = gdf[k].stool_age
		eegs = gdf[k].eeg_age

		y = subind[k.subject_id]
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

ylims!(ax_longsamples, -2, length(unique(long_sub.subject_id)) + 2)


# tweak some visuals

rowsize!(grid_longsamples, 3, Relative(5/6))
colsize!(figure1.layout, 1, Relative(1/3))
rowsize!(figure1.layout, 1, Relative(1/4))

hidedecorations!(ax_eeg_curves, ticks=false, ticklabels=false, label=false, minorticks=false)
hidedecorations!(ax_pcoa_spec, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_pcoa_func, ticks=false, ticklabels=false, label=false)
hidedecorations!(ax_longsamples, ticks=false, ticklabels=false, label=false)
hidexdecorations!.((ax_stool_hist, ax_eeg_hist))
hideydecorations!.((ax_stool_hist, ax_eeg_hist); ticks=false, ticklabels=false, label=false)
linkyaxes!(ax_eeg_hist, ax_stool_hist)
linkxaxes!(ax_eeg_hist, ax_stool_hist, ax_longsamples)

save("/home/kevin/Downloads/figure1-inprogress.png", figure1)
save("manuscript/mainfigures/figure1.svg", figure1)


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
				 yticks = (tickrange, [rich(g;color=colors_gstypes[gs_types_rev[g]]) for g in reverse(gssig)]),
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
				 yticks = (tickrange, [rich(g;color=colors_gstypes[gs_types_rev[g]]) for g in reverse(gssig)]),
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

save("/home/kevin/Downloads/figure2-inprogress.png", figure2)
save("manuscript/mainfigures/figure2.svg", figure2)

############
# Figure 3 #
############

figure3 = Figure(; size = (1050, 750))


grid_future_violins = GridLayout(figure3[1,1])
grid_futfsea_dots = GridLayout(figure3[1,2])
ax_future_violins = map(enumerate(ftps)) do (i, tp)
	ax = Axis(grid_future_violins[i,1]; ylabel = "age (months)", xticks = ([1,2], ["stool", "eeg"]),
		 xlabel=i == 3 ? "collection type" : "", title = "visit $(i == 3 ? "2" : "1") → visit $(i == 1 ? "2" : "3")")
	hidedecorations!(ax; label=false, ticks=false, ticklabels=false)
	ax
end

grid_futfsea_latency=GridLayout(grid_futfsea_dots[1,1])
ax_futlat = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_futfsea_latency[1,i];
				 xlabel = "z",
				 xticklabelsize= 10,
				 title = lab,
				 yticks = (tickrange, [rich(g;color=colors_gstypes[gs_types_rev[g]]) for g in reverse(gssig)]),
				 yticklabelsize=10,
				 )
		r = Axis(grid_futfsea_latency[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						   repeat(["v2→v3", "v1→v3", "v1→v2"]; outer=length(tickrange)))
				 )
		(l, r)
	end
end

grid_futfsea_amplitude=GridLayout(grid_futfsea_dots[1,2])
ax_futamp = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_futfsea_amplitude[1,i];
				 xlabel = "z",
				 xticklabelsize=10,
				 title = lab,
				 yticks = (tickrange, [rich(g;color=colors_gstypes[gs_types_rev[g]]) for g in reverse(gssig)]),
				 yticklabelsize=10,
				 )
		r = Axis(grid_futfsea_amplitude[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(["v2→v3", "v1→v3", "v1→v2"]; outer=length(tickrange)))
				 )
		(l, r)
	end
end

for axs in (ax_futlat, ax_futamp)
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
hideydecorations!(ax_futamp[1][1], ticks=false)

Legend(grid_future_violins[4,1], 
	   [[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 1:3],
		[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in reverse(5:7)]],
	   [["q < 0.01", "q < 0.1", "q < 0.2"], 
		reverse(["q < 0.2", "q < 0.1", "q < 0.01"])],
	   ["(-)", "(+)"];
	   orientation=:horizontal, tellheight=true, tellwidth=false, nbanks=3)

#-

for (i, df) in enumerate((v1v2, v1v3, v2v3))
	ax = ax_future_violins[i]
	clrs = [colors_timepoints[i == 3 ? 2 : 1][2], colors_timepoints[i == 1 ? 2 : 3][2]]
		    
	violin!(ax, repeat([1,2]; inner=size(df,1)), [df.stool_age; df.vep_age];
			color=repeat([(c, 0.3) for c in clrs]; inner=size(df, 1)))
	for row in eachrow(df)
		xs = [1,2] .+ rand(Normal(0,0.03))
		ys = [row.stool_age, row.vep_age]
		lines!(ax, xs, ys; color=:gray60, linestyle=:dash, linewidth=0.5)
		scatter!(ax, xs, ys; color=[(c, 0.3) for c in clrs], strokewidth=0.5, markersize=fsea_marker_size)
	end
end

#-


for (j, feat) in enumerate(filter(contains("latency"), eeg_features))
	gdf = groupby(futfsea_df, "eeg_feature")
    featstr = replace(feat, "peak_latency_"=>"", "_corrected"=>"c")
    @warn "$featstr"
    for (i, tp) in enumerate(ftps)
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
            violin!(ax_futlat[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_futlat[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_futlat[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
end


for (j, feat) in enumerate(filter(contains("amp"), eeg_features))
	gdf = groupby(futfsea_df, "eeg_feature")
    featstr = replace(feat, "peak_amp"=>"", "_corrected"=>"c")
    @warn "$featstr"
    for (i, tp) in enumerate(ftps)
        @info "$tp"

		df =  alllms[(; eeg_feature=feat, timepoint=tp)]
		allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        # lines!(ax_futlat, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
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
            violin!(ax_futamp[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_futamp[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_futamp[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
end

#-
colsize!(figure3.layout, 2, Relative(4/5))
linkyaxes!(ax_future_violins...)

save("/home/kevin/Downloads/figure3-inprogress.png", figure3)
save("manuscript/mainfigures/figure3.svg", figure3)

# ## Summary Ideas
#
# We really need some way to summarize results,
# for a possible figure 4.
# `./analyis/network.jl` has some code to generate
# edgelists that can be imported into cytoscape,
# which is one possibility. Here are some others.
#
# ### Volcano plots

volcano = Figure(; size=(700,900))

for (i, v) in enumerate(tps)
	df = subset(fsea_df, "timepoint"=> ByRow(==(v)))
	ax = Axis(volcano[i, 1]; xlabel="E.S.", ylabel="log(Q)")
	scatter!(df.es, -1 .* log.(df.q₀ .+ 0.0005); color= map(eachrow(df)) do row
		row.q₀ < 0.2 || return :gray
		colors_gstypes[gs_types_rev[row.geneset]]
	end)
	Label(volcano[i,0], v; tellheight=false)
end
	
Legend(volcano[1:3,2], [MarkerElement(; marker=:circle, color = colors_gstypes[t]) for t in keys(geneset_types)],
					   [String(t) for t in keys(geneset_types)]
)

save("/home/kevin/Downloads/volcano-fsea.png", current_figure())

#-

volcano = Figure(; size=(1500,800))

for (j, p) in enumerate(eeg_features)
	for (i, v) in enumerate(tps)
		df = subset(fsea_df, "timepoint"=> ByRow(==(v)), "eeg_feature"=> ByRow(==(p)))
		ax = Axis(volcano[i, j]; xlabel="E.S.", ylabel="log(Q)")
		scatter!(df.es, -1 .* log.(df.q₀ .+ 0.0005); color= map(eachrow(df)) do row
			row.q₀ < 0.2 || return :gray
			colors_gstypes[gs_types_rev[row.geneset]]
		end)
		Label(volcano[i,0], v; tellheight=false)
	end
	Label(volcano[0,j], p; tellwidth=false)
end
	
Legend(volcano[1:3,length(eeg_features)+1], [MarkerElement(; marker=:circle, color = colors_gstypes[t]) for t in keys(geneset_types)],
					   [String(t) for t in keys(geneset_types)]
)

save("/home/kevin/Downloads/volcano-fsea-peaks.png", current_figure())

# ### Vanja bars
#
# Another idea to summarize from Vanja, looking at counts of significant hits.

bars = Figure(; size=(900, 600))
axs = Axis[]

for (i, v) in enumerate(tps)
	df = subset(fsea_df, "timepoint"=> ByRow(==(v)))
	ax = Axis(bars[1,i]; ylabel="count", xticks=(1:length(eeg_features), eeg_features), xticklabelrotation=pi/2)
	toplot = mapreduce(vcat, enumerate(eeg_features)) do (j,f)
		 subdf = subset(df, "eeg_feature"=> ByRow(==(f)))
		 map(enumerate(collect(keys(geneset_types)))) do (k, gs)
			  n = count(<(0.2), subset(subdf, "geneset"=> ByRow(g-> g ∈ geneset_types[gs])).q₀)
			  (; feature=f, fidx = j, geneset=gs, gidx = k, n_sig=n)
		 end
	 end |> DataFrame

	barplot!(ax, toplot.fidx, toplot.n_sig; stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
	Label(bars[0, i], v; tellwidth=false)
	push!(axs, ax)

end

Legend(bars[1, length(tps)+1], [MarkerElement(; marker=:rect, color = colors_gstypes[t]) for t in keys(geneset_types)],
					   [String(t) for t in keys(geneset_types)]
)

linkyaxes!(axs...)

save("/home/kevin/Downloads/bars-fsea.png", current_figure())


#-

bars = Figure(; size=(900, 600))
axs = Axis[]

for (i, v) in enumerate(tps)
	df = subset(fsea_df, "timepoint"=> ByRow(==(v)))
	ax = Axis(bars[1,i]; ylabel="count", xticks=(1:length(eeg_features), eeg_features), xticklabelrotation=pi/2)
	toplot = mapreduce(vcat, enumerate(eeg_features)) do (j,f)
		 subdf = subset(df, "eeg_feature"=> ByRow(==(f)))
		 map(enumerate(collect(keys(geneset_types)))) do (k, gs)
			 npos = count(<(0.2), subset(subdf, "geneset"=> ByRow(g-> g ∈ geneset_types[gs]), "es"=> ByRow(>(0))).q₀)
			 nneg = count(<(0.2), subset(subdf, "geneset"=> ByRow(g-> g ∈ geneset_types[gs]), "es"=> ByRow(<(0))).q₀)
			 (; feature=f, fidx = j, geneset=gs, gidx = k, n_sig=(npos, nneg))
		 end
	end |> DataFrame

	barplot!(ax, toplot.fidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
	barplot!(ax, toplot.fidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
	Label(bars[0, i], v; tellwidth=false)

	push!(axs, ax)
end

Legend(bars[1, length(tps)+1], [MarkerElement(; marker=:rect, color = colors_gstypes[t]) for t in keys(geneset_types)],
					   [String(t) for t in keys(geneset_types)]
)

linkyaxes!(axs...)
save("/home/kevin/Downloads/bars-fsea-posneg.png", current_figure())
