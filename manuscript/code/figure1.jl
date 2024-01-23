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
	subset!(df, "feature"=> ByRow(f-> !contains(f, '|')))
	sample = replace(basename(f), r"(SEQ\d+).+"=> s"\1")
	df.sample .= sample
	df
    end
end

gfscomm = let 
    fs = unique(gfs.feature)
    ss = unique(gfs.sample)
    fsmap = Dict(f=> i for (i, f) in enumerate(fs))
    ssmap = Dict(s=> i for (i, s) in enumerate(ss))

    mat = spzeros(length(fs), length(ss))
    foreach(eachrow(gfs)) do row
	mat[fsmap[row.feature], ssmap[row.sample]] = row.abundance
    end
    CommunityProfile(mat, GeneFunction.(fs), MicrobiomeSample.(ss))
end


set!(gfscomm, mapreduce(get, vcat, (concurrent_3m, concurrent_6m, concurrent_12m)))
##

figure = Figure(; size=(800, 800))
ax = Axis(figure[1,1])
p = plot_pcoa!(ax, pcoa(gfscomm); color=get(gfscomm, :age))
Colorbar(figure[1,2], p; label="age (months)")

save("data/figures/pcoas/concurrent_all_functions_age.png", figure)

##

figure = Figure(; size=(1200,800))
ax1 = Axis(figure[1,1]; xticks=(1:2, ["stool", "eeg"]), ylabel="age (months)", title="3m stool -> 6m eeg")
ax2 = Axis(figure[1,2]; xticks=(1:2, ["stool", "eeg"]), ylabel="age (months)", title="3m stool -> 12m eeg")
ax3 = Axis(figure[1,3]; xticks=(1:2, ["stool", "eeg"]), ylabel="age (months)", title="3m stool -> 12m eeg")

let df = DataFrame(get(future_3m6m))
    xs1 = 1 .+ rand(Normal(0.0, 0.05), size(df, 1))
    xs2 = 2 .+ rand(Normal(0.0, 0.05), size(df, 1))
    violin!(ax1, repeat([1,2]; inner=size(df,1)), [df.age; df.eeg_age]; color=:gray80)
    scatter!(ax1, [xs1; xs2], [df.age; df.eeg_age]; color=(:gray50,0.5))
    foreach(i-> lines!(ax1, [xs1[i], xs2[i]], [df.age[i], df.eeg_age[i]]; color=:gray50), eachindex(xs1))
end
let df = DataFrame(get(future_3m12m))
    xs1 = 1 .+ rand(Normal(0.0, 0.05), size(df, 1))
    xs2 = 2 .+ rand(Normal(0.0, 0.05), size(df, 1))
    violin!(ax2, repeat([1,2]; inner=size(df,1)), [df.age; df.eeg_age]; color=:gray80)
    scatter!(ax2, [xs1; xs2], [df.age; df.eeg_age]; color=(:gray50,0.5))
    foreach(i-> lines!(ax2, [xs1[i], xs2[i]], [df.age[i], df.eeg_age[i]]; color=:gray50), eachindex(xs1))
end
let df = DataFrame(get(future_6m12m))
    xs1 = 1 .+ rand(Normal(0.0, 0.05), size(df, 1))
    xs2 = 2 .+ rand(Normal(0.0, 0.05), size(df, 1))
    violin!(ax3, repeat([1,2]; inner=size(df,1)), [df.age; df.eeg_age]; color=:gray80)
    scatter!(ax3, [xs1; xs2], [df.age; df.eeg_age]; color=(:gray50,0.5))
    foreach(i-> lines!(ax3, [xs1[i], xs2[i]], [df.age[i], df.eeg_age[i]]; color=:gray50), eachindex(xs1))
end
linkyaxes!(ax1, ax2, ax3)
save("data/figures/future_age_dists.png", figure)

##

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
    on="subject"), "subject", r"3m", r"6m", r"12m")

tps = ("3m", "6m", "12m")
ldf = DataFrame()
for row in eachrow(wide_sub)
    for tp in tps
	stool_age = row["seqprep_$tp"]
	eeg_age = row["eeg_$tp"]
	push!(ldf, (; subject=row.subject, timepoint=tp, stool_age, eeg_age); cols=:union)
    end
end
subset!(ldf, AsTable(["stool_age", "eeg_age"])=> ByRow(nt-> !all(ismissing, nt)))
transform!(ldf,AsTable(["stool_age", "eeg_age"])=> ByRow(nt-> minimum(skipmissing(values(nt))))=> "minage") 
sort!(ldf, "minage")
subind = Dict(s=> i for (i,s) in enumerate(unique(ldf.subject)))

gdf = groupby(ldf, "subject")

figure = Figure(; size=(500, 1200))
ax = Axis(figure[1,1]; xlabel="age (months)", ylabel = "subject")
stoolc = :firebrick 
eegc = :dodgerblue
bothc = :slateblue

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
	ismissing(s) && return eegc
	ismissing(e) && return stoolc
	return bothc
    end
    scatter!(ax, xs, fill(y, length(xs)); color=cs)
 
 
    lines!(ax, [extrema(skipmissing([stools;eegs]))...], [y,y]; linestyle=:dash, color = :gray)
    for s in filter(!ismissing, stools)
	lines!(ax, [s,s], [y+0.4, y-0.4]; color=stoolc)
    end
    for e in filter(!ismissing, eegs)
	lines!(ax, [e,e], [y+0.4, y-0.4]; color=eegc)
    end
end
save("data/figures/subject_longitudinal.png", figure)
