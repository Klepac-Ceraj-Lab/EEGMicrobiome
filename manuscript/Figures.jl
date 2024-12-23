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

using VKCComputing
using EEGMicrobiome
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


tps = ("v1", "v2", "v3")
ftps = ("v1v2", "v1v3", "v2v3")

mdata = load_cohorts()

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



##

eeg_features = "peak_latency_" .* ["N1", "P1_corrected", "N2_corrected"]
eeg_features = [eeg_features; replace.(eeg_features, "latency" => "amp")]

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
  CommunityProfile(abundances(mpa), features(mpa), MicrobiomeSample.(replace.(samplenames(mpa), r"_S\d++_profile" => "")))
end

set!(taxprofiles, select(mdata, "seqprep" => "sample", Cols(:)))

unirefs_by_sample = let
  files = String.(skipmissing(mdata.genefamilies))
  mapreduce(vcat, files) do f
    df = CSV.read(f, DataFrame)
    rename!(df, ["feature", "abundance"])
    subset!(df, "feature" => ByRow(f -> !contains(f, '|')))
    sample = replace(basename(f), r"(SEQ\d+).+" => s"\1")
    df.sample .= sample
    df
  end
end

unirefs = let
  fs = unique(unirefs_by_sample.feature)
  ss = unique(unirefs_by_sample.sample)
  fsmap = Dict(f => i for (i, f) in enumerate(fs))
  ssmap = Dict(s => i for (i, s) in enumerate(ss))

  mat = spzeros(length(fs), length(ss))
  foreach(eachrow(unirefs_by_sample)) do row
    mat[fsmap[row.feature], ssmap[row.sample]] = row.abundance
  end
  CommunityProfile(mat, GeneFunction.(fs), MicrobiomeSample.(ss))
end


set!(unirefs, select(mdata, "seqprep" => "sample", Cols(:)))

unirefs_pco = pcoa(unirefs)
species_pco = pcoa(taxprofiles)

fsea_df = CSV.read("data/outputs/fsea/concurrent_consolidated_fsea.csv", DataFrame)
futfsea_df = CSV.read("data/outputs/fsea/future_consolidated_fsea.csv", DataFrame)

geneset_types = (;
  neurotransmitters=[
    "GABA synthesis",
    "GABA degradation",
    "Glutamate synthesis",
    "Glutamate degradation",
  ],
  nt_metabolism=[
    # tryptophan/serotonin stuff
    "Tryptophan synthesis",
    "Tryptophan degradation",
    "Quinolinic acid synthesis",
    "Quinolinic acid degradation",
    # dopamine stuff
    "DOPAC synthesis",
    "DOPAC degradation",
  ],
  scfa=[
    "Acetate synthesis",
    "Acetate degradation",
    "Propionate synthesis",
    "Propionate degradation",
    "Butyrate synthesis",
    "Butyrate degradation",
    "Isovaleric acid synthesis",
    "Isovaleric acid degradation",
  ],
  other=[
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
gs_types_rev = Dict(gs => t for t in keys(geneset_types) for gs in geneset_types[t])

geneset_order = vcat(geneset_types...)
gssig = intersect(geneset_order,
  unique(subset(vcat(fsea_df, futfsea_df), "q₀" => ByRow(<(0.2))).geneset)
)
gsidx = geneset_index(gssig, geneset_order)

alllms = groupby(mapreduce(vcat, Iterators.product(eeg_features, [tps..., ftps...])) do (feat, tp)
    df = CSV.read("data/outputs/lms/$(feat)_$(tp)_lms.csv", DataFrame)
    rename!(df, "Name" => "eeg_feature")
    df.timepoint .= tp
    df
  end, ["eeg_feature", "timepoint"])

for gs in gssig
  @info gs
  transform!(alllms, "feature" => ByRow(u -> replace(u, "UniRef90_" => "") ∈ na_map[gs]) => gs; ungroup=false)
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
# ### Theming
#
# First, we want to set the size of physical units
# (Makie uses typical CSS conventions, see https://docs.makie.org/v0.21/how-to/match-figure-size-font-sizes-and-dpi)/

inch = 96
pt = 4/3
cm = inch / 2.54

# Then, we set the typical size of text in different figure contexts

mytheme = Theme(;
    font = "Liberation Sans",
    Axis = (; 
        xticklabelsize = 8pt,
        xlabelsize = 10pt,
        xgridvisible = false,
        yticklabelsize = 8pt,
        ylabelsize = 10pt,
        ygridvisible = false,
        titlesize = 10pt
    ),
    Label = (;
        fontsize = 10pt
    ),
    Legend = (;
        labelsize = 8pt,
        titlesize = 10pt
    )
)

update_theme!(mytheme)

#
# ### Colors

colormap_age = :viridis
colors_timepoints = [tp => c for (tp, c) in zip(tps, cgrad(colormap_age)[[0.0, 0.4, 0.8]])]
colors_sampletype = [st => c for (st, c) in zip(["stool", "eeg", "both"], cgrad(:tab10; categorical=true)[[1, 4, 5]])]
colors_eeg_kind = [kd => c for (kd, c) in zip(["amp", "latency"], cgrad(:tab10; categorical=true)[[2, 3]])]
colors_eeg_peaks = [pk => c for (pk, c) in zip(["N1", "P1c", "N2c"], cgrad(:Accent_4; categorical=true)[2:4])]

colormap_sig = :PuOr
colors_sig = cgrad(colormap_sig, 11; rev=true, categorical=true)[[1, 3, 4, 8, 9, 11]]
insert!(colors_sig, 4, colorant"white")

colors_gstypes = NamedTuple(gst => c for (gst, c) in zip(keys(geneset_types), cgrad(:tab10; categorical=true)[[10, 3, 2, 9]]))

# ### Tables
#
# #### Table 2

table2 = subset(fsea_df, "timepoint" => ByRow(t -> t ∈ tps),
  "q₀" => ByRow(<(0.2)),
  "geneset" => ByRow(g -> g ∈ geneset_order)
)
transform!(table2, "eeg_feature" => ByRow(f -> begin
  (_, p, n...) = split(f, "_")
  n = first(n)
  return (; feature=p, peak=n)
end) => ["feature", "peak"]
)
transform!(table2, "timepoint" => ByRow(x -> Dict("v1" => 1, "v2" => 2, "v3" => 3)[x]) => "visit")
select!(table2, "visit", "geneset", "feature", "peak", "es" => "E.S.", "q₀" => "q value")
sort!(table2, [
  "visit",
  order("geneset"; by=x -> findfirst(==(x), geneset_order)),
  "feature",
  # order("peak"; by = x-> findfirst(==(x), ("N1", "P1", "N2")))
  "E.S."
]
)

for v in 1:3
  CSV.write("/home/kevin/Downloads/Table2_v$v.csv", subset(table2, "visit" => ByRow(==(v))))
end


# #### Table 3

table3 = subset(futfsea_df,
  "timepoint" => ByRow(t -> t ∈ ftps),
  "q₀" => ByRow(<(0.2)),
  "geneset" => ByRow(g -> g ∈ geneset_order)
)
transform!(table3, "eeg_feature" => ByRow(f -> begin
  (_, p, n...) = split(f, "_")
  n = first(n)
  return (; feature=p, peak=n)
end) => ["feature", "peak"]
)
transform!(table3, "timepoint" => ByRow(x -> begin
  (st, vep) = match(r"^v(\d+)v(\d+)$", x).captures
  return (; stool=parse(Int, st), VEP=parse(Int, vep))
end) => ["stool", "VEP"]
)
select!(table3, "stool", "VEP", "geneset", "feature", "peak", "es" => "E.S.", "q₀" => "q value")
sort!(table3, [
  "stool",
  "VEP",
  order("geneset"; by=x -> findfirst(==(x), geneset_order)),
  "feature",
  # order("peak"; by = x-> findfirst(==(x), ("N1", "P1", "N2")))
  "E.S."
]
)
for (s, v) in zip((1, 1, 2), (2, 3, 3))
  CSV.write("/home/kevin/Downloads/Table3_$(s)_$(v).csv", subset(table3, "stool" => ByRow(==(s)), "VEP" => ByRow(==(v))))
end

# ### Main figures
#
# #### Figure 1

figure1 = Figure(; size=(10inch, 7.5inch))

grid_cohort = GridLayout(figure1[1, 1:2])
grid_longsamples = GridLayout(figure1[2:3, 1])
grid_eeg_curves = GridLayout(figure1[3, 2])

grid_pcoas = GridLayout(figure1[2, 2])

ax_cohort = Axis(grid_cohort[1, 1]; aspect=DataAspect(), alignmode=Outside())
# Legend(grid_cohort[2,1],
# 	[MarkerElement(; marker=:rect, color=c[2]) for c in colors_timepoints],
# 	["visit 1", "visit 2", "visit 3"];
# 	orientation=:horizontal, tellwidth=false, tellheight=true
# )

ax_eeg_curves = Axis(grid_eeg_curves[1, 1]; xlabel="time (ms) relative to stimulus onset",
  ylabel="voltage (μV)",
  yminorticksvisible=true,
  yminorticks=IntervalsBetween(5)
)
ax_pcoa_spec = Axis(grid_pcoas[1, 1])
ax_pcoa_func = Axis(grid_pcoas[1, 2])

ax_longsamples = Axis(grid_longsamples[3, 1];
  xlabel="age (months)",
  xticks=0:3:18,
  limits=((0.0, 19.0), nothing),
  ylabel="subject"
)
ax_eeg_hist = Axis(grid_longsamples[1, 1]; ylabel="density", title="eeg")
ax_stool_hist = Axis(grid_longsamples[2, 1]; ylabel="density", title="stool")


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

let datage = data(subset(mdata, "eeg_age" => ByRow(!ismissing))) * mapping(:eeg_age => "age (months)", color="visit")
    draw!(ax_eeg_hist, datage * AlgebraOfGraphics.density(), scales(Color= (; palette=colors_timepoints)))
end

let datage = data(subset(mdata, "stool_age" => ByRow(!ismissing))) * mapping(:stool_age => "age (months)", color="visit")
    draw!(ax_stool_hist, datage * AlgebraOfGraphics.density(), scales(Color = (; palette = colors_timepoints)))
end

mpl = draw!(ax_eeg_curves, datmean * visual(Lines), scales(Color=(; palette = colors_timepoints)))
dpl = draw!(ax_eeg_curves, datlower * visual(Lines; linestyle=:dash), scales(Color=(; palette = colors_timepoints)))
draw!(ax_eeg_curves, datupper * visual(Lines; linestyle=:dash), scales(Color=(; palette = colors_timepoints)))

Legend(grid_eeg_curves[1, 1],
  [LineElement(; color=:gray), LineElement(; color=:gray, linestyle=:dash)],
  ["mean", "+/- S.E."];
  tellwidth=false, halign=:right, valign=:bottom, margin=(10, 10, 10, 10)
)

Legend(grid_eeg_curves[1, 1],
  [LineElement(; color=c[2]) for c in colors_timepoints],
  ["Visit $i" for i in 1:3];
  tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10)
)
# ##### PCoAs


plt_pcoa_spec = plot_pcoa!(ax_pcoa_spec, species_pco; color=get(taxprofiles, :stool_age), colormap=colormap_age)
plt_pcoa_func = plot_pcoa!(ax_pcoa_func, unirefs_pco; color=get(unirefs, :stool_age), colormap=colormap_age)

Colorbar(grid_pcoas[1, 3];
  limits=extrema(skipmissing(mdata.stool_age)), label="Age (months)", colormap=colormap_age,
)

# ##### Longitudinal


subind = Dict(s => i for (i, s) in enumerate(unique(long_sub.subject_id)))

samples_colors = [s => c for (s, c) in zip(
  ("stool", "eeg", "both"),
  (colorant"dodgerblue", colorant"firebrick", colorant"slateblue"))
]

let gdf = groupby(long_sub, "subject_id")
  for k in keys(gdf)
    stools = gdf[k].stool_age
    eegs = gdf[k].eeg_age

    y = subind[k.subject_id]
    xs = map(zip(stools, eegs)) do (s, e)
      ismissing(s) && return e
      ismissing(e) && return s
      return mean([e, s])
    end
    cs = map(zip(stools, eegs)) do (s, e)
      ismissing(s) && return colors_sampletype[2][2]
      ismissing(e) && return colors_sampletype[1][2]
      return colors_sampletype[3][2]
    end
    scatter!(ax_longsamples, xs, fill(y, length(xs)); color=cs)


    lines!(ax_longsamples, [extrema(skipmissing([stools; eegs]))...], [y, y]; linestyle=:dash, color=:gray)
    for s in filter(!ismissing, stools)
      lines!(ax_longsamples, [s, s], [y + 0.4, y - 0.4]; color=colors_sampletype[1][2])
    end
    for e in filter(!ismissing, eegs)
      lines!(ax_longsamples, [e, e], [y + 0.4, y - 0.4]; color=colors_sampletype[2][2])
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

rowsize!(grid_longsamples, 3, Relative(5 / 6))
colsize!(figure1.layout, 1, Relative(1 / 3))
rowsize!(figure1.layout, 1, Relative(1 / 4))

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
figure1

# ##### Figure 2
#
# ###### Setup

function topspec(profile, n=10)
    specsum = [sum(abundances(profile[i, :])) for i in 1:size(profile, 1)]
    rank = invperm(sortperm(specsum; rev=true))
    topn = findall(r-> r ∈ 1:n, rank)
    other = sum(abundances(profile)[Not(topn), :]; dims=1)
    abmat = vcat(abundances(profile)[topn, :], other)
    CommunityProfile(abmat, [features(profile)[topn]; [Microbiome.Taxon("other")]], samples(profile))
end

v1_top = topspec(taxprofiles[:, [s for s in v1.seqprep]], 8)
v2_top = topspec(taxprofiles[:, [s for s in v2.seqprep]], 8)
v3_top = topspec(taxprofiles[:, [s for s in v3.seqprep]], 8)

v1_top_dm = Microbiome.braycurtis(v1_top)
v2_top_dm = Microbiome.braycurtis(v2_top)
v3_top_dm = Microbiome.braycurtis(v3_top)

v1_top_hcl = hclust(v1_top_dm; linkage=:average, branchorder=:optimal)
v2_top_hcl = hclust(v2_top_dm; linkage=:average, branchorder=:optimal)
v3_top_hcl = hclust(v3_top_dm; linkage=:average, branchorder=:optimal)

colors_species = [sp => c for (sp, c) in zip(
    [filter(!=("other"), union(featurenames.([v1_top, v2_top, v3_top])...)); ["other"]],
    cgrad(:Paired_12; categorical=true)[[[1:10...]; [12, 11]]])
]

# ###### Plotting

figure2 = Figure(; size=(10inch,7.5inch));

top_grid = GridLayout(figure2[1,1])
bottom_grid = GridLayout(figure2[2,1])
eeg_scatters = GridLayout(top_grid[1,1])
tax_abundances = GridLayout(bottom_grid[1,1])
fsea_grid = GridLayout(top_grid[1, 2])
fsea_violins = GridLayout(bottom_grid[1,2])


for (j,feat) in enumerate(filter(f-> contains(f, "amp"), eeg_features))
    ax = Axis(eeg_scatters[j, 2]; xticklabelsvisible = j == 3)
    for (i, tp) in enumerate(tps)
        df = subset(mdata, "cohort_$tp" => identity)
        scatter!(ax, df[!, "eeg_age"], df[!, feat]; color = colors_timepoints[i][2])
    end
end

for (j,feat) in enumerate(filter(f-> contains(f, "lat"), eeg_features))
    ax = Axis(eeg_scatters[j, 4]; xticklabelsvisible = j == 3)
    for (i, tp) in enumerate(tps)
        df = subset(mdata, "cohort_$tp" => identity)
        scatter!(ax, df[!, "eeg_age"], df[!, feat]; color = colors_timepoints[i][2])
    end
end


for (j, feat) in enumerate(filter(f-> contains(f, "amp"), eeg_features))
    peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
    Label(eeg_scatters[j,0], peak; tellwidth=true, tellheight=false)
end
Label(eeg_scatters[1:3,1], "voltage (μV)"; rotation=π/2, tellwidth=true, tellheight=false)
Label(eeg_scatters[1:3,3], "latency(ms)"; rotation=π/2, tellwidth=true, tellheight=false)
Label(eeg_scatters[4, 2:4], "age (months)"; tellheight=true, tellwidth=false)

v1_ab_ax = Axis(tax_abundances[1,1]; xticksvisible=false, xticklabelsvisible=false)
v2_ab_ax = Axis(tax_abundances[1,2]; xticksvisible=false, xticklabelsvisible=false)
v3_ab_ax = Axis(tax_abundances[1,3]; xticksvisible=false, xticklabelsvisible=false)

let prof = v1_top[:, v1_top_hcl.order]
    df = DataFrame(species = [s for s in repeat(featurenames(prof), size(prof, 2))],
                   sample = [s for s in repeat(samplenames(prof), inner=size(prof, 1))],
                   abundance = [a for a in vec(abundances(prof))]
    )
    spec = data(df) * mapping(:sample => presorted, :abundance; stack = :species, color=:species)
    draw!(v1_ab_ax, spec * visual(BarPlot), scales(Color = (; palette = colors_species)))
end

let prof = v2_top[:, v2_top_hcl.order]
    df = DataFrame(species = [s for s in repeat(featurenames(prof), size(prof, 2))],
                   sample = [s for s in repeat(samplenames(prof), inner=size(prof, 1))],
                   abundance = [a for a in vec(abundances(prof))]
    )
    spec = data(df) * mapping(:sample => presorted, :abundance; stack = :species, color=:species)
    draw!(v2_ab_ax, spec * visual(BarPlot), scales(Color = (; palette = colors_species)))
end

let prof = v3_top[:, v3_top_hcl.order]
    df = DataFrame(species = [s for s in repeat(featurenames(prof), size(prof, 2))],
                   sample = [s for s in repeat(samplenames(prof), inner=size(prof, 1))],
                   abundance = [a for a in vec(abundances(prof))]
    )
    spec = data(df) * mapping(:sample => presorted, :abundance; stack = :species, color=:species)
    draw!(v3_ab_ax, spec * visual(BarPlot), scales(Color = (; palette = colors_species)))
end

for (i, v) in enumerate(["v1", "v2", "v3"])
    Label(tax_abundances[2,i], v; tellheight=true, tellwidth=false)
end

Label(tax_abundances[3, 1:3], "sample"; tellheight=true, tellwidth=false)
Label(tax_abundances[1, 0], "abundance (%)"; rotation= π/2, tellheight=false, tellwidth=true)

tightlimits!.([v1_ab_ax, v2_ab_ax, v3_ab_ax])

Legend(tax_abundances[4,:],
       [MarkerElement(; marker=:rect, color=c[2]) for c in colors_species],
       map(colors_species) do c
            c[1] == "other" && return rich(c[1])
            rich(replace(c[1], "_"=> " "), font=:italic)
       end;
       rowgap=-3, colgap=1, nbanks=5, orientation=:horizontal,
       padding=(3f0, 3f0, 3f0, 3f0))

# rowsize!(figure2.layout, 2, Relative(3/7))

#-

let fsea_bars_axs = Axis[]

    for (j, feat) in enumerate(filter(f-> contains(f, "amp"),  eeg_features))
        df = subset(fsea_df, "eeg_feature" => ByRow(==(feat)))
        ax = Axis(fsea_grid[j, 2]; xticks=(1:length(tps), [tps...]))
        toplot = mapreduce(vcat, enumerate(tps)) do (i,tp)
            subdf = subset(df, "timepoint" => ByRow(==(tp)))
            map(enumerate(collect(keys(geneset_types)))) do (k, gs)
                npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
                nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
                (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
            end
        end |> DataFrame

        peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
        
        barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        hlines!(ax, [0.]; color=:black, linewidth=0.5 )

        push!(fsea_bars_axs, ax)
    end

    linkyaxes!(fsea_bars_axs...)
end

let fsea_bars_axs = Axis[]

    for (j, feat) in enumerate(filter(f-> contains(f, "lat"),  eeg_features))
        df = subset(fsea_df, "eeg_feature" => ByRow(==(feat)))
        ax = Axis(fsea_grid[j, 4]; xticks=(1:length(tps), [tps...]))
        toplot = mapreduce(vcat, enumerate(tps)) do (i,tp)
            subdf = subset(df, "timepoint" => ByRow(==(tp)))
            map(enumerate(collect(keys(geneset_types)))) do (k, gs)
                npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
                nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
                (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
            end
        end |> DataFrame

        barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        hlines!(ax, [0.]; color=:black, linewidth=0.5 )


        push!(fsea_bars_axs, ax)
    end

    linkyaxes!(fsea_bars_axs...)
end

Label(fsea_grid[1:3,1], "Significant genesets w/amplitude (N)"; rotation=π/2, tellwidth=true, tellheight=false)
Label(fsea_grid[1:3,3], "Significant genesets w/latency (N)"; rotation=π/2, tellwidth=true, tellheight=false)

Legend(fsea_grid[4, :],
    [MarkerElement(; marker=:rect, color=colors_gstypes[t]) for t in keys(geneset_types)],
    ["Neurotransmitters", "Amino acid metabolism", "SCFAs", "other"];
    orientation=:horizontal, padding=(3f0, 3f0, 3f0, 3f0), nbanks=2
)


fsea_violin1 = Axis(fsea_violins[1,1]; ylabel="E.S.", xticks=(1:3, ["v1", "v2", "v3"]))
fsea_enrichment1 = GridLayout(fsea_violins[2,1])
fsea_violin2 = Axis(fsea_violins[1,2]; ylabel="E.S.", xticks=(1:3, ["v1", "v2", "v3"]))
fsea_enrichment2 = GridLayout(fsea_violins[2,2])

gs_interval = 6
tp_interval = 0.5
dotsxlim = 3
fsea_marker_size = 6

let (j, feat, gs) = (1, "peak_latency_N2_corrected", "Glutamate synthesis")
    gdf = groupby(fsea_df, "eeg_feature")
    featstr = replace(feat, "peak_latency_" => "", "_corrected" => "c")
    @warn "$featstr"
    for (i, tp) in enumerate(tps)
        @info "$tp"

        df = alllms[(; eeg_feature=feat, timepoint=tp)]
        allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        row = first(subset(subdf, "geneset"=> ByRow(==(gs))))
        yidx = df[!, row.geneset]
        xpos = i
        ys = df.z[yidx]
        xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
        ymed = median(ys)

        cidx = row.q₀ > 0.2 ? 4 :       # not significant
                row.q₀ < 0.01 ? 1 :       # quite significant
                row.q₀ < 0.1 ? 2 : 3 # somewhat significant / not very significant
        c = ymed < allymed ? colors_sig[cidx] : colors_sig[8-cidx]
        violin!(fsea_violin1, fill(xpos, length(ys)), ys; width=tp_interval * 1.5, color=(c, 0.4))
        scatter!(fsea_violin1, xs, ys; color=c, strokewidth=0.5, markersize=fsea_marker_size)
        lines!(fsea_violin1, [xpos - tp_interval / 2, xpos + tp_interval / 2], [ymed, ymed]; color=c)

        fr = GridLayout(fsea_enrichment1[i, 1])
        fsea_res = FeatureSetEnrichments.FSEAResult(row.q₀, row.nfeatures, eval(Meta.parse(row.ranks)))
        plot!(fr, fsea_res; linecolor = c == colors_sig[4] ? colorant"lightgray" : c, xticks = WilkinsonTicks(3;k_max=4), ylabelvisible=false)

    end
end

let (j, feat, gs) = (1, "peak_amp_P1_corrected", "Menaquinone synthesis")
    gdf = groupby(fsea_df, "eeg_feature")
    featstr = replace(feat, "peak_latency_" => "", "_corrected" => "c")
    @warn "$featstr"
    for (i, tp) in enumerate(tps)
        @info "$tp"

        df = alllms[(; eeg_feature=feat, timepoint=tp)]
        allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        row = first(subset(subdf, "geneset"=> ByRow(==(gs))))
        yidx = df[!, row.geneset]
        xpos = i
        ys = df.z[yidx]
        xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
        ymed = median(ys)

        cidx = row.q₀ > 0.2 ? 4 :       # not significant
                row.q₀ < 0.01 ? 1 :       # quite significant
                row.q₀ < 0.1 ? 2 : 3 # somewhat significant / not very significant
        c = ymed < allymed ? colors_sig[cidx] : colors_sig[8-cidx]
        violin!(fsea_violin2, fill(xpos, length(ys)), ys; width=tp_interval * 1.5, color=(c, 0.4))
        scatter!(fsea_violin2, xs, ys; color=c, strokewidth=0.5, markersize=fsea_marker_size)
        lines!(fsea_violin2, [xpos - tp_interval / 2, xpos + tp_interval / 2], [ymed, ymed]; color=c)

        fr = GridLayout(fsea_enrichment2[i, 1])
        fsea_res = FeatureSetEnrichments.FSEAResult(row.q₀, row.nfeatures, eval(Meta.parse(row.ranks)))
        (gr, ax1, ax2) = plot!(fr, fsea_res; linecolor = c == colors_sig[4] ? colorant"lightgray" : c, xticks = WilkinsonTicks(3; k_max=4), ylabelvisible=false)
    end
end

colsize!(top_grid, 1, Relative(2/3))
colsize!(bottom_grid, 1, Relative(2/3))


save("/home/kevin/Downloads/figure2-inprogress.png", figure2)
save("manuscript/mainfigures/figure2.svg", figure2)
figure2

# ##### Figure S2

figureS2 = Figure(; size=(1100, 600))


grid_fsea_dots = GridLayout(figureS2[1, 1])
legend_block = GridLayout(figureS2[2, 1])
# grid_fsea_heatmaps = GridLayout(figureS2[1,2])
# grid_bug_heatmaps = GridLayout(figureS2[2, 1])

gs_interval = 6
tp_interval = 1.5
dotsxlim = 3
fsea_marker_size = 6

grid_fsea_latency = GridLayout(grid_fsea_dots[1, 1])
ax_lat = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
  map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
    l = Axis(grid_fsea_latency[1, i];
      xlabel="z",
      title=lab,
      yticks=(tickrange, [rich("$g ", rich("■"; color=colors_gstypes[gs_types_rev[g]])) for g in reverse(gssig)]),
    )
    r = Axis(grid_fsea_latency[1, i];
      yaxisposition=:right,
      yticklabelsize=6,
      yticks=(sort([collect(tickrange);
          collect(tickrange) .- tp_interval;
          collect(tickrange) .+ tp_interval]),
        repeat(string.(3:-1:1); outer=length(tickrange)))
    )
    (l, r)
  end
end

grid_fsea_amplitude = GridLayout(grid_fsea_dots[1, 2])
ax_amp = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
  map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
    l = Axis(grid_fsea_amplitude[1, i];
      xlabel="z",
      title=lab,
      yticks=(tickrange, [rich("$g ", rich("■"; color=colors_gstypes[gs_types_rev[g]])) for g in reverse(gssig)]),
    )
    r = Axis(grid_fsea_amplitude[1, i];
      yaxisposition=:right,
      yticklabelsize=6,
      yticks=(sort([collect(tickrange);
          collect(tickrange) .- tp_interval;
          collect(tickrange) .+ tp_interval]),
        repeat(string.(3:-1:1); outer=length(tickrange)))
    )
    (l, r)
  end
end

for axs in (ax_lat, ax_amp)
  for (i, (axl, axr)) in enumerate(axs)
    ylims!.((axl, axr), gs_interval - 2 * tp_interval, gs_interval * length(gssig) + 2 * tp_interval)
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

Legend(legend_block[1, 1],
  [MarkerElement(; marker=:rect, color=colors_gstypes[g]) for g in keys(geneset_types)][[1, 3, 2, 4]],
  ["Neurotransmitters", "Amino acid metabolism", "SCFAs", "other"][[1, 3, 2, 4]];
  tellheight=true, tellwidth=false, nbanks=2
)


Legend(legend_block[1, 2],
  [[MarkerElement(; marker=:circle, color=colors_sig[i]) for i in 1:3],
    [MarkerElement(; marker=:circle, color=colors_sig[i]) for i in 5:7]],
  [["q < 0.01", "q < 0.1", "q < 0.2"],
    ["q < 0.2", "q < 0.1", "q < 0.01"]],
  ["(-)", "(+)"];
  orientation=:horizontal, tellheight=true, tellwidth=false
)

#-


for (j, feat) in enumerate(filter(contains("latency"), eeg_features))
  gdf = groupby(fsea_df, "eeg_feature")
  featstr = replace(feat, "peak_latency_" => "", "_corrected" => "c")
  @warn "$featstr"
  for (i, tp) in enumerate(tps)
    @info "$tp"

    df = alllms[(; eeg_feature=feat, timepoint=tp)]
    allymed = median(df.z)

    subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
    for row in eachrow(subdf)
      yidx = df[!, row.geneset]
      xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
      ys = df.z[yidx]
      xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
      ymed = median(ys)

      cidx = row.q₀ > 0.2 ? 4 :       # not significant
             row.q₀ < 0.01 ? 1 :       # quite significant
             row.q₀ < 0.1 ? 2 : 3 # somewhat significant / not very significant
      c = ymed < allymed ? colors_sig[cidx] : colors_sig[8-cidx]
      violin!(ax_lat[j][1], fill(xpos, length(ys)), ys; width=tp_interval * 1.5, color=(c, 0.4), orientation=:horizontal)
      scatter!(ax_lat[j][1], ys, xs; color=c, strokewidth=0.5, markersize=fsea_marker_size)
      lines!(ax_lat[j][1], [ymed, ymed], [xpos - tp_interval / 2, xpos + tp_interval / 2]; color=c)
    end
  end
end


for (j, feat) in enumerate(filter(contains("amp"), eeg_features))
  gdf = groupby(fsea_df, "eeg_feature")
  featstr = replace(feat, "peak_amp" => "", "_corrected" => "c")
  @warn "$featstr"
  for (i, tp) in enumerate(tps)
    @info "$tp"

    df = alllms[(; eeg_feature=feat, timepoint=tp)]
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

      cidx = row.q₀ > 0.2 ? 4 :       # not significant
             row.q₀ < 0.01 ? 1 :       # quite significant
             row.q₀ < 0.1 ? 2 : 3 # somewhat significant / not very significant
      c = ymed < allymed ? colors_sig[cidx] : colors_sig[8-cidx]
      violin!(ax_amp[j][1], fill(xpos, length(ys)), ys; width=tp_interval * 1.5, color=(c, 0.4), orientation=:horizontal)
      scatter!(ax_amp[j][1], ys, xs; color=c, strokewidth=0.5, markersize=fsea_marker_size)
      lines!(ax_amp[j][1], [ymed, ymed], [xpos - tp_interval / 2, xpos + tp_interval / 2]; color=c)
    end
  end
end

colsize!(legend_block, 1, Relative(1 / 8))

save("/home/kevin/Downloads/figureS2-inprogress.png", figureS2)
save("manuscript/mainfigures/figureS2.svg", figureS2)

############
# Figure 3 #
############

figure3 = Figure(; size=(6inch,4inch));

grid_future_violins = GridLayout(figure3[1, 1]; alignmode=Outside())
grid_fsea_bars = GridLayout(figure3[1,2])

ax_future_violins = map(enumerate(ftps)) do (i, tp)
  ax = Axis(grid_future_violins[i, 1]; ylabel="age (months)", xticks=([1, 2], ["stool", "eeg"]),
    xlabel=i == 3 ? "collection type" : "", title="visit $(i == 3 ? "2" : "1") → visit $(i == 1 ? "2" : "3")")
end

linkyaxes!(ax_future_violins...)

for (i, df) in enumerate((v1v2, v1v3, v2v3))
  ax = ax_future_violins[i]
  clrs = [colors_timepoints[i == 3 ? 2 : 1][2], colors_timepoints[i == 1 ? 2 : 3][2]]

  violin!(ax, repeat([1, 2]; inner=size(df, 1)), [df.stool_age; df.vep_age];
    color=repeat([(c, 0.3) for c in clrs]; inner=size(df, 1)))
  for row in eachrow(df)
    xs = [1, 2] .+ rand(Normal(0, 0.03))
    ys = [row.stool_age, row.vep_age]
    lines!(ax, xs, ys; color=:gray60, linestyle=:dash, linewidth=0.5)
    scatter!(ax, xs, ys; color=[(c, 0.3) for c in clrs], strokewidth=0.5, markersize=fsea_marker_size)
  end
end

let fsea_bars_axs = Axis[]

    for (j, feat) in enumerate(filter(f-> contains(f, "amp"),  eeg_features))
        df = subset(futfsea_df, "eeg_feature" => ByRow(==(feat)))
        ax = Axis(grid_fsea_bars[j, 2]; xticks=(1:length(ftps), [ftps...]))
        toplot = mapreduce(vcat, enumerate(ftps)) do (i,tp)
            subdf = subset(df, "timepoint" => ByRow(==(tp)))
            map(enumerate(collect(keys(geneset_types)))) do (k, gs)
                npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
                nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
                (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
            end
        end |> DataFrame

        peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
        
        barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        # Label(grid_fsea_bars[0, j], peak; tellwidth=false)

        hlines!(ax, [0.]; color=:black, linewidth=0.5 )
        push!(fsea_bars_axs, ax)
    end

    linkyaxes!(fsea_bars_axs...)
end

let fsea_bars_axs = Axis[]

    for (j, feat) in enumerate(filter(f-> contains(f, "lat"),  eeg_features))
        df = subset(futfsea_df, "eeg_feature" => ByRow(==(feat)))
        ax = Axis(grid_fsea_bars[j, 4]; xticks=(1:length(ftps), [ftps...]))
        toplot = mapreduce(vcat, enumerate(ftps)) do (i,tp)
            subdf = subset(df, "timepoint" => ByRow(==(tp)))
            map(enumerate(collect(keys(geneset_types)))) do (k, gs)
                npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
                nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
                (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
            end
        end |> DataFrame

        barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
        hlines!(ax, [0.]; color=:black, linewidth=0.5 )

        push!(fsea_bars_axs, ax)
    end

    linkyaxes!(fsea_bars_axs...)
end

Label(grid_fsea_bars[1:3,1], "Significant genesets w/amplitude (N)"; rotation=π/2, tellwidth=true, tellheight=false)
Label(grid_fsea_bars[1:3,3], "Significant genesets w/latency (N)"; rotation=π/2, tellwidth=true, tellheight=false)

for (j, feat) in enumerate(filter(f-> contains(f, "amp"), eeg_features))
    peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
    Label(grid_fsea_bars[j,0], peak; tellwidth=true, tellheight=false)
end

Legend(grid_fsea_bars[4, :],
    [MarkerElement(; marker=:rect, color=colors_gstypes[t]) for t in keys(geneset_types)],
    ["Neurotransmitters", "Amino acid metabolism", "SCFAs", "other"];
    orientation=:horizontal, padding=(3f0, 3f0, 3f0, 3f0), nbanks=2
)

colsize!(figure3.layout, 1, Relative(1/3))

save("/home/kevin/Downloads/figure3-inprogress.png", figure3)
save("manuscript/mainfigures/figure3.svg", figure3)

figure3
#-

figureS3 = Figure(; size=(1050, 750))


grid_future_violins = GridLayout(figureS3[1:2, 1]; alignmode=Outside())
grid_futfsea_dots = GridLayout(figureS3[1, 2]; alignmode=Outside())
legend_block = GridLayout(figureS3[2, 2])

ax_future_violins = map(enumerate(ftps)) do (i, tp)
  ax = Axis(grid_future_violins[i, 1]; ylabel="age (months)", xticks=([1, 2], ["stool", "eeg"]),
    xlabel=i == 3 ? "collection type" : "", title="visit $(i == 3 ? "2" : "1") → visit $(i == 1 ? "2" : "3")")
  hidedecorations!(ax; label=false, ticks=false, ticklabels=false)
  ax
end

grid_futfsea_latency = GridLayout(grid_futfsea_dots[1, 1])
ax_futlat = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
  map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
    l = Axis(grid_futfsea_latency[1, i];
      xlabel="z",
      xticklabelsize=10,
      title=lab,
      yticks=(tickrange, [rich("$g ", rich("■"; color=colors_gstypes[gs_types_rev[g]])) for g in reverse(gssig)]),
      yticklabelsize=10,
    )
    r = Axis(grid_futfsea_latency[1, i];
      yaxisposition=:right,
      yticklabelsize=6,
      yticks=(sort([collect(tickrange);
          collect(tickrange) .- tp_interval;
          collect(tickrange) .+ tp_interval]),
        repeat(["v2→v3", "v1→v3", "v1→v2"]; outer=length(tickrange)))
    )
    (l, r)
  end
end

grid_futfsea_amplitude = GridLayout(grid_futfsea_dots[1, 2])
ax_futamp = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
  map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
    l = Axis(grid_futfsea_amplitude[1, i];
      xlabel="z",
      xticklabelsize=10,
      title=lab,
      yticks=(tickrange, [rich(g; color=colors_gstypes[gs_types_rev[g]]) for g in reverse(gssig)]),
      yticklabelsize=10,
    )
    r = Axis(grid_futfsea_amplitude[1, i];
      yaxisposition=:right,
      yticklabelsize=6,
      yticks=(sort([collect(tickrange);
          collect(tickrange) .- tp_interval;
          collect(tickrange) .+ tp_interval]),
        repeat(["v2→v3", "v1→v3", "v1→v2"]; outer=length(tickrange)))
    )
    (l, r)
  end
end

for axs in (ax_futlat, ax_futamp)
  for (i, (axl, axr)) in enumerate(axs)
    ylims!.((axl, axr), gs_interval - 2 * tp_interval, gs_interval * length(gssig) + 2 * tp_interval)
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


#-

for (i, df) in enumerate((v1v2, v1v3, v2v3))
  ax = ax_future_violins[i]
  clrs = [colors_timepoints[i == 3 ? 2 : 1][2], colors_timepoints[i == 1 ? 2 : 3][2]]

  violin!(ax, repeat([1, 2]; inner=size(df, 1)), [df.stool_age; df.vep_age];
    color=repeat([(c, 0.3) for c in clrs]; inner=size(df, 1)))
  for row in eachrow(df)
    xs = [1, 2] .+ rand(Normal(0, 0.03))
    ys = [row.stool_age, row.vep_age]
    lines!(ax, xs, ys; color=:gray60, linestyle=:dash, linewidth=0.5)
    scatter!(ax, xs, ys; color=[(c, 0.3) for c in clrs], strokewidth=0.5, markersize=fsea_marker_size)
  end
end

#-


for (j, feat) in enumerate(filter(contains("latency"), eeg_features))
  gdf = groupby(futfsea_df, "eeg_feature")
  featstr = replace(feat, "peak_latency_" => "", "_corrected" => "c")
  @warn "$featstr"
  for (i, tp) in enumerate(ftps)
    @info "$tp"

    df = alllms[(; eeg_feature=feat, timepoint=tp)]
    allymed = median(df.z)

    subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
    for row in eachrow(subdf)
      yidx = df[!, row.geneset]
      xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
      ys = df.z[yidx]
      xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
      ymed = median(ys)

      cidx = row.q₀ > 0.2 ? 4 :       # not significant
             row.q₀ < 0.01 ? 1 :       # quite significant
             row.q₀ < 0.1 ? 2 : 3 # somewhat significant / not very significant
      c = ymed < allymed ? colors_sig[cidx] : colors_sig[8-cidx]
      violin!(ax_futlat[j][1], fill(xpos, length(ys)), ys; width=tp_interval * 1.5, color=(c, 0.4), orientation=:horizontal)
      scatter!(ax_futlat[j][1], ys, xs; color=c, strokewidth=0.5, markersize=fsea_marker_size)
      lines!(ax_futlat[j][1], [ymed, ymed], [xpos - tp_interval / 2, xpos + tp_interval / 2]; color=c)
    end
  end
end


for (j, feat) in enumerate(filter(contains("amp"), eeg_features))
  gdf = groupby(futfsea_df, "eeg_feature")
  featstr = replace(feat, "peak_amp" => "", "_corrected" => "c")
  @warn "$featstr"
  for (i, tp) in enumerate(ftps)
    @info "$tp"

    df = alllms[(; eeg_feature=feat, timepoint=tp)]
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

      cidx = row.q₀ > 0.2 ? 4 :       # not significant
             row.q₀ < 0.01 ? 1 :       # quite significant
             row.q₀ < 0.1 ? 2 : 3 # somewhat significant / not very significant
      c = ymed < allymed ? colors_sig[cidx] : colors_sig[8-cidx]
      violin!(ax_futamp[j][1], fill(xpos, length(ys)), ys; width=tp_interval * 1.5, color=(c, 0.4), orientation=:horizontal)
      scatter!(ax_futamp[j][1], ys, xs; color=c, strokewidth=0.5, markersize=fsea_marker_size)
      lines!(ax_futamp[j][1], [ymed, ymed], [xpos - tp_interval / 2, xpos + tp_interval / 2]; color=c)
    end
  end
end

#-
Legend(legend_block[1, 1],
  [MarkerElement(; marker=:rect, color=colors_gstypes[g]) for g in keys(geneset_types)][[1, 3, 2, 4]],
  ["Neurotransmitters", "Amino acid metabolism", "SCFAs", "other"][[1, 3, 2, 4]];
  tellheight=true, tellwidth=false, nbanks=2
)


Legend(legend_block[1, 2],
  [[MarkerElement(; marker=:circle, color=colors_sig[i]) for i in 1:3],
    [MarkerElement(; marker=:circle, color=colors_sig[i]) for i in 5:7]],
  [["q < 0.01", "q < 0.1", "q < 0.2"],
    ["q < 0.2", "q < 0.1", "q < 0.01"]],
  ["(-)", "(+)"];
  orientation=:horizontal, tellheight=true, tellwidth=false
)

colsize!(figureS3.layout, 2, Relative(4 / 5))
linkyaxes!(ax_future_violins...)

save("/home/kevin/Downloads/figureS3-inprogress.png", figureS3)
save("manuscript/mainfigures/figureS3.svg", figureS3)

# ## Summary Ideas
#
# We really need some way to summarize results,
# for a possible figure 4.
# `./analyis/network.jl` has some code to generate
# edgelists that can be imported into cytoscape,
# which is one possibility. Here are some others.
#
# ### Volcano plots

volcano = Figure(; size=(700, 900))

for (i, v) in enumerate(tps)
  df = subset(fsea_df, "timepoint" => ByRow(==(v)))
  ax = Axis(volcano[i, 1]; xlabel="E.S.", ylabel="log(Q)")
  scatter!(df.es, -1 .* log.(df.q₀ .+ 0.0005); color=map(eachrow(df)) do row
    row.q₀ < 0.2 || return :gray
    colors_gstypes[gs_types_rev[row.geneset]]
  end)
  Label(volcano[i, 0], v; tellheight=false)
end

Legend(volcano[1:3, 2], [MarkerElement(; marker=:circle, color=colors_gstypes[t]) for t in keys(geneset_types)],
  [String(t) for t in keys(geneset_types)]
)

save("/home/kevin/Downloads/volcano-fsea.png", current_figure())

#-

volcano = Figure(; size=(1500, 800))

for (j, p) in enumerate(eeg_features)
  for (i, v) in enumerate(tps)
    df = subset(fsea_df, "timepoint" => ByRow(==(v)), "eeg_feature" => ByRow(==(p)))
    ax = Axis(volcano[i, j]; xlabel="E.S.", ylabel="log(Q)")
    scatter!(df.es, -1 .* log.(df.q₀ .+ 0.0005); color=map(eachrow(df)) do row
      row.q₀ < 0.2 || return :gray
      colors_gstypes[gs_types_rev[row.geneset]]
    end)
    Label(volcano[i, 0], v; tellheight=false)
  end
  Label(volcano[0, j], p; tellwidth=false)
end

Legend(volcano[1:3, length(eeg_features)+1], [MarkerElement(; marker=:circle, color=colors_gstypes[t]) for t in keys(geneset_types)],
  [String(t) for t in keys(geneset_types)]
)

save("/home/kevin/Downloads/volcano-fsea-peaks.png", current_figure())

# ### Vanja bars
#
# Another idea to summarize from Vanja, looking at counts of significant hits.

bars = Figure(; size=(900, 600));
axs = Axis[]

for (j, feat) in enumerate(filter(f-> contains(f, "amp"),  eeg_features))
    df = subset(fsea_df, "eeg_feature" => ByRow(==(feat)))
    ax = Axis(bars[1, j]; ylabel="count", xticks=(1:length(tps), [tps...]))
    toplot = mapreduce(vcat, enumerate(tps)) do (i,tp)
        subdf = subset(df, "timepoint" => ByRow(==(tp)))
        map(enumerate(collect(keys(geneset_types)))) do (k, gs)
            npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
            nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
            (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
        end
    end |> DataFrame

    peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
    
    barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    Label(bars[0, j], peak; tellwidth=false)

    push!(axs, ax)
end

linkyaxes!(axs...)

axs = Axis[]

for (j, feat) in enumerate(filter(f-> contains(f, "lat"),  eeg_features))
    df = subset(fsea_df, "eeg_feature" => ByRow(==(feat)))
    ax = Axis(bars[2, j]; ylabel="count", xticks=(1:length(tps), [tps...]))
    toplot = mapreduce(vcat, enumerate(tps)) do (i,tp)
        subdf = subset(df, "timepoint" => ByRow(==(tp)))
        map(enumerate(collect(keys(geneset_types)))) do (k, gs)
            npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
            nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
            (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
        end
    end |> DataFrame

    barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])

    push!(axs, ax)
end

linkyaxes!(axs...)

Label(bars[1, 0], "Amplitude"; tellwidth=true, tellheight=false)
Label(bars[2, 0], "Latency"; tellwidth=true, tellheight=false)

Legend(bars[:, 4], [MarkerElement(; marker=:rect, color=colors_gstypes[t]) for t in keys(geneset_types)],
  [String(t) for t in keys(geneset_types)]
)
bars

save("/home/kevin/Downloads/bars-fsea-by-peak.png", bars)

#-

bars = Figure(; size=(900, 600));
axs = Axis[]

for (j, feat) in enumerate(filter(f-> contains(f, "amp"),  eeg_features))
    df = subset(futfsea_df, "eeg_feature" => ByRow(==(feat)))
    ax = Axis(bars[1, j]; ylabel="count", xticks=(1:length(ftps), [ftps...]))
    toplot = mapreduce(vcat, enumerate(ftps)) do (i,tp)
        subdf = subset(df, "timepoint" => ByRow(==(tp)))
        map(enumerate(collect(keys(geneset_types)))) do (k, gs)
            npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
            nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
            (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
        end
    end |> DataFrame

    peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
    
    barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    Label(bars[0, j], peak; tellwidth=false)

    push!(axs, ax)
end

linkyaxes!(axs...)

axs = Axis[]

for (j, feat) in enumerate(filter(f-> contains(f, "lat"),  eeg_features))
    df = subset(futfsea_df, "eeg_feature" => ByRow(==(feat)))
    ax = Axis(bars[2, j]; ylabel="count", xticks=(1:length(ftps), [ftps...]))
    toplot = mapreduce(vcat, enumerate(ftps)) do (i,tp)
        subdf = subset(df, "timepoint" => ByRow(==(tp)))
        map(enumerate(collect(keys(geneset_types)))) do (k, gs)
            npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
            nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
            (; feature=feat, tidx=i, geneset=gs, gidx=k, n_sig=(npos, nneg))
        end
    end |> DataFrame

    barplot!(ax, toplot.tidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    barplot!(ax, toplot.tidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])

    push!(axs, ax)
end

linkyaxes!(axs...)

Label(bars[1, 0], "Amplitude"; tellwidth=true, tellheight=false)
Label(bars[2, 0], "Latency"; tellwidth=true, tellheight=false)

Legend(bars[:, 4], [MarkerElement(; marker=:rect, color=colors_gstypes[t]) for t in keys(geneset_types)],
  [String(t) for t in keys(geneset_types)]
)
bars

save("/home/kevin/Downloads/bars-futfsea-by-peak.png", bars)
#-

bars = Figure(; size=(900, 600))
axs = Axis[]

for (i, v) in enumerate(tps)
    df = subset(fsea_df, "timepoint" => ByRow(==(v)))
    ax = Axis(bars[1, i]; ylabel="count", xticks=(1:length(eeg_features), eeg_features), xticklabelrotation=pi / 2)
    toplot = mapreduce(vcat, enumerate(eeg_features)) do (j, f)
        subdf = subset(df, "eeg_feature" => ByRow(==(f)))
        map(enumerate(collect(keys(geneset_types)))) do (k, gs)
        npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
        nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
        (; feature=f, fidx=j, geneset=gs, gidx=k, n_sig=(npos, nneg))
        end
    end |> DataFrame

    barplot!(ax, toplot.fidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    barplot!(ax, toplot.fidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
    Label(bars[0, i], v; tellwidth=false)

    push!(axs, ax)
end

Legend(bars[1, length(tps)+1], [MarkerElement(; marker=:rect, color=colors_gstypes[t]) for t in keys(geneset_types)],
  [String(t) for t in keys(geneset_types)]
)

linkyaxes!(axs...)
save("/home/kevin/Downloads/bars-fsea-posneg.png", current_figure())

#-

bars = Figure(; size=(900, 600))

ampaxs = Axis[]

for (i, v) in enumerate(ftps)
  df = subset(futfsea_df, "timepoint" => ByRow(==(v)))
  ax = Axis(bars[1, i]; ylabel="count", xticks=(1:length(eeg_features), eeg_features), xticklabelrotation=pi / 2)
    toplot = mapreduce(vcat, enumerate(filter(f-> contains(f, "amp", eeg_features)))) do (j, f)
    subdf = subset(df, "eeg_feature" => ByRow(==(f)))
    map(enumerate(collect(keys(geneset_types)))) do (k, gs)
      npos = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(>(0))).q₀)
      nneg = count(<(0.2), subset(subdf, "geneset" => ByRow(g -> g ∈ geneset_types[gs]), "es" => ByRow(<(0))).q₀)
      (; feature=f, fidx=j, geneset=gs, gidx=k, n_sig=(npos, nneg))
    end
  end |> DataFrame

  barplot!(ax, toplot.fidx, getindex.(toplot.n_sig, 1); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
  barplot!(ax, toplot.fidx, -1 .* getindex.(toplot.n_sig, 2); stack=toplot.gidx, color=[colors_gstypes[gs] for gs in toplot.geneset])
  Label(bars[0, i], v; tellwidth=false)

  push!(ampaxs, ax)
end

Legend(bars[1, length(tps)+1], [MarkerElement(; marker=:rect, color=colors_gstypes[t]) for t in keys(geneset_types)],
  [String(t) for t in keys(geneset_types)]
)

linkyaxes!(ampaxs...)
save("/home/kevin/Downloads/bars-fsea-future-posneg.png", current_figure())

# ### Age / EEG scatters

amp_scatters = Figure(; size=(900,900), title="Amplitudes");

for (i, tp) in enumerate(tps), (j,feat) in enumerate(filter(f-> contains(f, "amp"), eeg_features))
    ax = Axis(amp_scatters[j, i]; xlabel="age (months)", ylabel="voltage (μV)")
    df = subset(mdata, "cohort_$tp" => identity)
    scatter!(ax, df[!, "eeg_age"], df[!, feat])
end

for (i, tp) in enumerate(tps)
    Label(amp_scatters[0,i], tp; tellwidth=false)
end

for (j, feat) in enumerate(filter(f-> contains(f, "amp"), eeg_features))
    peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
    Label(amp_scatters[j,0], peak; tellwidth=true, tellheight=false)
end

amp_scatters

save("/home/kevin/Downloads/eeg_amp_scatters_concurrent.png", amp_scatters)

lat_scatters = Figure(; size=(900,900), title="Latencies");

for (i, tp) in enumerate(tps), (j,feat) in enumerate(filter(f-> contains(f, "lat"), eeg_features))
    ax = Axis(lat_scatters[j, i]; xlabel="age (months)", ylabel="latency (ms)")
    df = subset(mdata, "cohort_$tp" => identity)
    scatter!(ax, df[!, "eeg_age"], df[!, feat])
end

for (i, tp) in enumerate(tps)
    Label(lat_scatters[0,i], tp; tellwidth=false)
end

for (j, feat) in enumerate(filter(f-> contains(f, "lat"), eeg_features))
    peak = replace(feat, r"peak_latency_([NP12]+)(_corrected)?" => s"\1")
    Label(lat_scatters[j,0], peak; tellwidth=true, tellheight=false)
end

lat_scatters

save("/home/kevin/Downloads/eeg_lat_scatters_concurrent.png", lat_scatters)

#-

amp_scatters = Figure(; size=(900,900), title="Amplitudes");

for (j,feat) in enumerate(filter(f-> contains(f, "amp"), eeg_features))
    ax = Axis(amp_scatters[j, 1]; xlabel="age (months)", ylabel="voltage (μV)")
    for (i, tp) in enumerate(tps)
        df = subset(mdata, "cohort_$tp" => identity)
        scatter!(ax, df[!, "eeg_age"], df[!, feat]; color = colors_timepoints[i][2])
    end
end

for (j, feat) in enumerate(filter(f-> contains(f, "amp"), eeg_features))
    peak = replace(feat, r"peak_amp_([NP12]+)(_corrected)?" => s"\1")
    Label(amp_scatters[j,0], peak; tellwidth=true, tellheight=false)
end

amp_scatters

save("/home/kevin/Downloads/eeg_amp_scatters_tp_combined.png", amp_scatters)

lat_scatters = Figure(; size=(900,900), title="Latencies");

for (j,feat) in enumerate(filter(f-> contains(f, "lat"), eeg_features))
    ax = Axis(lat_scatters[j, 1]; xlabel="age (months)", ylabel="latency (ms)")
    for (i, tp) in enumerate(tps)
        df = subset(mdata, "cohort_$tp" => identity)
        scatter!(ax, df[!, "eeg_age"], df[!, feat]; color = colors_timepoints[i][2])
    end
end

for (j, feat) in enumerate(filter(f-> contains(f, "lat"), eeg_features))
    peak = replace(feat, r"peak_latency_([NP12]+)(_corrected)?" => s"\1")
    Label(lat_scatters[j,0], peak; tellwidth=true, tellheight=false)
end

lat_scatters

save("/home/kevin/Downloads/eeg_lat_scatters_tp_combined.png", lat_scatters)

# ### Ordinations

summary_ord = Figure(; size=(900, 699))
ax = Axis(summary_ord[1, 1])



######################
# Manuscript numbers #
######################


let
  lead = "Stool samples and EEG were collected at up to 3 visits in the first 18 months of life (Figure 1A, B, Table 1;"
  v1all = subset(mdata, AsTable(["cohort_v1", "cohort_v1v2_stool", "cohort_v1v3_stool"]) => ByRow(any))
  v2all = subset(mdata, AsTable(["cohort_v2", "cohort_v1v2_vep", "cohort_v2v3_stool"]) => ByRow(any))
  v3all = subset(mdata, AsTable(["cohort_v3", "cohort_v1v3_vep", "cohort_v2v3_vep"]) => ByRow(any))

  v1ages = transform(v1all, AsTable(["stool_age", "eeg_age"]) => ByRow(nt -> mean(skipmissing(values(nt)))) => "mean_age").mean_age
  v2ages = transform(v2all, AsTable(["stool_age", "eeg_age"]) => ByRow(nt -> mean(skipmissing(values(nt)))) => "mean_age").mean_age
  v3ages = transform(v3all, AsTable(["stool_age", "eeg_age"]) => ByRow(nt -> mean(skipmissing(values(nt)))) => "mean_age").mean_age

  println(lead, " visit 1, N = $(length(v1ages)), age $(round(mean(v1ages); digits=1)) ± $(round(std(v1ages); digits=1)) months,",
    " visit 2, N = $(length(v2ages)), age $(round(mean(v2ages); digits=1)) ± $(round(std(v2ages); digits=1)) months,",
    " visit 3, N = $(length(v3ages)), age $(round(mean(v3ages); digits=1)) ± $(round(std(v3ages); digits=1)) months)"
  )
end

println("(A) Study design; participants (N=", length(unique(mdata.subject_id)), ")")

println("Of the ", count(k -> length(na_map[k]) > 5, keys(na_map)), " genesets tested, ",
  length(unique(fsea_df.geneset)), " had sufficient representation to test, and of those, ",
  length(unique(subset(fsea_df, "q₀" => ByRow(<(0.2))).geneset)), " were significantly associated with at least one")

println("To investigate this, we performed FSEA on stool samples collected at visit 1 with VEP measured at visit 2 ",
  "(age at stool collection = ", round(mean(subset(mdata, "cohort_v1v2_stool").stool_age); digits=1), " ± ",
  round(std(subset(mdata, "cohort_v1v2_stool").stool_age); digits=1), " months, ",
  "age at VEP = ", round(mean(subset(mdata, "cohort_v1v2_vep").eeg_age); digits=1), " ± ",
  round(std(subset(mdata, "cohort_v1v2_vep").eeg_age); digits=1), " months)"
)

println("or visit 3 ",
  "(age at stool collection = ", round(mean(subset(mdata, "cohort_v1v3_stool").stool_age); digits=1), " ± ",
  round(std(subset(mdata, "cohort_v1v3_stool").stool_age); digits=1), " months, ",
  "age at VEP = ", round(mean(subset(mdata, "cohort_v1v3_vep").eeg_age); digits=1), " ± ",
  round(std(subset(mdata, "cohort_v1v3_vep").eeg_age); digits=1), " months)"
)


println("as well as visit 2 stool samples with visit 3 VEP ",
  "(age at stool collection = ", round(mean(subset(mdata, "cohort_v2v3_stool").stool_age); digits=1), " ± ",
  round(std(subset(mdata, "cohort_v2v3_stool").stool_age); digits=1), " months, ",
  "age at VEP = ", round(mean(subset(mdata, "cohort_v2v3_vep").eeg_age); digits=1), " ± ",
  round(std(subset(mdata, "cohort_v2v3_vep").eeg_age); digits=1), " months)"
)

println("primarily at the 3rd visit (mean age ", round(subset(mdata, "cohort_v3" => identity).stool_age |> mean, digits=1), " months)")
println("also associated with the latency of the P1 peak at visit 1 (mean age ", round(subset(mdata, "cohort_v1" => identity).stool_age |> mean, digits=1), " months)")

using MultivariateStats

cor(EEGMicrobiome.loadings(species_pco, 1), mdata.stool_age)
cor(EEGMicrobiome.loadings(unirefs_pco, 1), mdata.stool_age)

