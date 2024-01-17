using EEGMicrobiome
using VKCComputing
using CSV
using DataFrames
using CairoMakie
using Microbiome
using BiobakeryUtils

concurrent_3m = load_cohort("concurrent_3m")
concurrent_6m = load_cohort("concurrent_6m")
concurrent_12m = load_cohort("concurrent_12m")
future_3m6m = load_cohort("future_3m6m")
future_3m12m = load_cohort("future_3m12m")
future_6m12m = load_cohort("future_6m12m")


