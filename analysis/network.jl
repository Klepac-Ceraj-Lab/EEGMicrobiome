# # Network analysis
#
# Building network diagram with neuroactive gene, VEP features, and species

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
using Combinatorics
using StatsBase
using Distances
using Base.Threads

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


concurrent_species = commjoin(concurrent_3m, concurrent_6m, concurrent_12m)[:, samplenames(concurrent_unirefs)]
concurrent_species = filter(t-> taxrank(t) === :species, concurrent_species)
concurrent_species = commjoin([concurrent_species[:, [s]] for s in samplenames(concurrent_unirefs)]...)

@assert samplenames(concurrent_unirefs) == samplenames(concurrent_species)

nadict = FeatureSetEnrichments.get_neuroactive_unirefs()
na_reverse = Dict("UniRef90_$f" => k for k in keys(nadict) for f in nadict[k])
nacomm = let fs = Set("UniRef90_$u" for u in reduce(union, values(nadict)))
    filter(f-> name(f) in fs, concurrent_unirefs)
end
nacomm = filter(f-> mean(>(0), abundances(nacomm[f, :])) > 0.05, nacomm)

mdata = DataFrame(get(concurrent_species))
eeg_featues = ("peak_amp_N1", "peak_amp_P1_corrected", "peak_amp_N2_corrected")
eeg_featuers = vcat(eeg_features, replace.(eeg_features, "amp"=>"latency"))

allfeatures = vcat(featurenames(concurrent_species), featurenames(nacomm), eeg_features)
assmat = zeros(length(allfeatures), length(allfeatures))
spl = nfeatures(concurrent_species)


@threads for comb in Combinatorics.combinations(collect(eachindex(allfeatures)), 2) |> collect
    (i,j) = comb
    i < j || continue

    assmat[i, j] = cor(
        abundances(i <= spl ? concurrent_species[i, :] : nacomm[i-spl, :]),
        abundances(j <= spl ? concurrent_species[j, :] : nacomm[j-spl, :]);
        dims = 2
    ) |> only 
end


#-

using Graphs
using MetaGraphsNext
using GraphIO
using GraphMakie

assgraph = MetaGraph(
    Graph();
    label_type=String,
    vertex_data_type=Dict{String,String},
    edge_data_type=Float64,
    graph_data="Microbiome association graph"
)
for f in allfeatures
    assgraph[f] = Dict()
end

for i in eachindex(allfeatures), j in eachindex(allfeatures)
    i < j || continue
    iname = i <= spl ? featurenames(concurrent_species)[i] : featurenames(nacomm)[i-spl]
    assgraph[iname]["node_type"] = contains(iname, "UniRef") ? "function" : "taxon"

    jname = j <= spl ? featurenames(concurrent_species)[j] : featurenames(nacomm)[j-spl]

    assgraph[iname, jname] = assmat[i,j]
end

edges_df = DataFrame(map(edges(assgraph)) do edge
    i = edge.src
    j = edge.dst
    iname = i <= spl ? featurenames(concurrent_species)[i] : featurenames(nacomm)[i-spl]
    jname = j <= spl ? featurenames(concurrent_species)[j] : featurenames(nacomm)[j-spl]
    (; source_node=iname, dest_node=jname, i, j, edge_weight=assgraph[iname,jname])
    
end)

transform!(edges_df, "edge_weight"=> ByRow(w-> w < 0. ? "negative" : "positive") => "edgetype")
transform!(edges_df, AsTable(["source_node", "dest_node", "edgetype"]) => ByRow(nt-> string(nt.source_node, " (", nt.edgetype, ") ", nt.dest_node))=> "edge_name")

nodetable = DataFrame(name = allfeatures,
                      kind = map(f-> contains(f, "UniRef") ? "function" : "species", allfeatures),
                      prevalence = map(i -> only(prevalence(i <= spl ? concurrent_species[i,:] : nacomm[i - spl, :])), eachindex(allfeatures)),
                      neuroactive = map(f-> get(na_reverse, f, "n/a"), allfeatures)
)
