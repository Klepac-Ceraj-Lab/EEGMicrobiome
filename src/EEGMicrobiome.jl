module EEGMicrobiome

export load_eeg,
       load_microbiome,
       load_taxonomic_profiles!,
       load_functional_profiles!

using CSV
using DataFrames
using Chain
using VKCComputing
using Preferences
using BiobakeryUtils: metaphlan_profiles
using ThreadsX
using HypothesisTests: MannWhitneyUTest, pvalue
using MultipleTesting: adjust, BenjaminiHochberg
using GLM

include("data_loading.jl")
include("lms.jl")
include("bugs.jl")

end