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
using XLSX
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
using AlgebraOfGraphics
using SparseArrays

# Next, we'll load the data.
# The `load_cohort()` fucntion is written in `src/data_loading.jl`.

concurrent_3m = load_cohort("concurrent_3m")
concurrent_6m = load_cohort("concurrent_6m")
concurrent_12m = load_cohort("concurrent_12m")
future_3m6m = load_cohort("future_3m6m")
future_3m12m = load_cohort("future_3m12m")
future_6m12m = load_cohort("future_6m12m")

# Data for the EEG timeseries panel in Figure 1
# comes from separate analysis by Emma Margolis.

timeseries = mapreduce(vcat, ("3m", "6m", "12m")) do tp
    df = DataFrame(XLSX.readtable("data/alltimepoints_timeseries_figure.xlsx", "$tp Timeseries Calculation"; infer_eltypes=true))[:, 1:6]
    rename!(df, ["ms", "mean", "se", "lower", "upper", "std"])
    df.timepoint .= tp
    df[end, 1] = df[end-1, 1] + 1
    dropmissing!(df)
end

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

figure1 = Figure(; size = (1200,1200))

grid_cohort = GridLayout(figure1[1,1])
grid_eeg_curves = GridLayout(figure1[2,1])
grid_pcoas = GridLayout(figure1[1:2, 2])
grid_longsamples = GridLayout(figure1[1:2, 3])

#- 
ax_cohort = Axis(grid_cohort[1,1]; aspect = DataAspect(), alignmode=Outside())
hidedecorations!(ax_cohort)
hidespines!(ax_cohort)

image!(ax_cohort, rotr90(load("manuscript/mainfigures/eegfig1a.png")))



## 

datmean = data(timeseries) * mapping(:ms, :mean, color=:timepoint)
datlower = data(timeseries) * mapping(:ms, :lower, color=:timepoint)
datupper = data(timeseries) * mapping(:ms, :upper, color=:timepoint)

cs = [tp=> c for (tp, c) in zip(("3m", "6m", "12m"), cgrad(:batlow, 3; categorical=true))]

sub = Axis(grid_eeg_curves[1,1]; xlabel="time (ms) relative to stimulus onset",
		     ylabel="voltage (μV)")

mpl = draw!(sub, datmean * visual(Lines); palettes=(; color=cs))
dpl = draw!(sub, datlower * visual(Lines; linestyle=:dash); palettes=(; color=cs))
draw!(sub, datupper * visual(Lines; linestyle=:dash); palettes=(; color=cs))

Legend(grid_eeg_curves[1,2], [[LineElement(; color=cs[i][2]) for i in 1:3], [LineElement(; color=:gray), LineElement(; color=:gray, linestyle=:dash)]],
		  [["visit 1", "visit 2", "visit 3"], ["mean", "+/- S.E."]], ["visit", "value"]
) 

