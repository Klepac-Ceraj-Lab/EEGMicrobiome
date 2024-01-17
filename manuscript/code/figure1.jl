using EEGMicrobiome
using VKCComputing
using CSV
using DataFrames
using CairoMakie
using Microbiome
using BiobakeryUtils
using CairoMakie

concurrent_3m = load_cohort("concurrent_3m")
concurrent_6m = load_cohort("concurrent_6m")
concurrent_12m = load_cohort("concurrent_12m")
future_3m6m = load_cohort("future_3m6m")
future_3m12m = load_cohort("future_3m12m")
future_6m12m = load_cohort("future_6m12m")


## 

figure = Figure(; size=(800,800))
ax = Axis(figure[1,1])
let comm = commjoin(concurrent_3m, concurrent_6m, concurrent_12m)
    p = plot_pcoa!(ax, pcoa(comm); color=get(comm, :age))
    Colorbar(figure[1,2], p; label = "age (months)")
end

save("data/figures/pcoas/concurrent_all_age.png", figure)


gfs = let ss = mapreduce(samplenames, vcat, (concurrent_3m, concurrent_6m, concurrent_12m))
    files = filter(f-> contains(basename(f), "genefamilies.tsv") && replace(basename(f), r"(SEQ\d+).+"=> s"\1") ∈ ss,
		readdir("/grace/sequencing/processed/mgx/humann/main/"; join=true)
    )
    mapreduce(vcat, files) do f
	df = CSV.read(f, DataFrame)
	rename!(df, ["feature", "abundance"])
	sample = replace(basename(f), r"(SEQ\d+).+"=> s"\1")
	df.sample .= sample
    end
end


