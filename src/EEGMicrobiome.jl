module EEGMicrobiome

export load_cohorts,
       get_cohort

export plot_pcoa!,
       topbugs,
       grouptop,
       plotbugs!
    
export geneset_index,
       format_eeg_str

using CSV
using DataFrames
using Chain
using VKCComputing
using Preferences
using BiobakeryUtils: metaphlan_profiles
using Microbiome
using ThreadsX
using HypothesisTests: MannWhitneyUTest, pvalue
using MultipleTesting: adjust, BenjaminiHochberg
using GLM
using MultivariateStats
using CairoMakie

include("data_loading.jl")
include("lms.jl")
include("bugs.jl")
include("plotting.jl")

end
