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



@assert samplenames(taxprofiles) == samplenames(unirefs)
v1_idx = get(taxprofiles, :cohort_v1)
v2_idx = get(taxprofiles, :cohort_v2)
v3_idx = get(taxprofiles, :cohort_v3)
v1v2_idx = get(taxprofiles, :cohort_v1v2_stool)
v1v3_idx = get(taxprofiles, :cohort_v1v3_stool)
v2v3_idx = get(taxprofiles, :cohort_v2v3_stool)

using PERMANOVA

function permanovas(comm::CommunityProfile, metadatums; n = 1000, mdlabels = String.(metadatums))
    permdf = mapreduce(vcat, zip(metadatums, mdlabels)) do (md, lab)    
        com_md = get(comm, md)
        hasmd = findall(!ismissing, com_md)
        df = DataFrame(test = com_md[hasmd])
        disallowmissing!(df)
        size(df, 1) < 20 && @warn "Very small number of not missing $lab"

        p = permanova(df, abundances(comm)[:, hasmd]', BrayCurtis, @formula(1 ~ test), n)

        return (; metadatum = lab, varexpl=varexpl(p), pvalue=pvalue(p))
    end

    return DataFrame(permdf)
end

varexpl(p::PERMANOVA.PSummary) = p.results[1, 3] * 100
FeatureSetEnrichments.pvalue(p::PERMANOVA.PSummary) = p.results[1, 5]


pnovs = DataFrame()
for (idx, label) in zip((v1_idx, v2_idx, v3_idx, v1v2_idx, v1v3_idx, v2v3_idx),
                        ("v1_idx", "v2_idx", "v3_idx", "v1v2_idx", "v1v3_idx", "v2v3_idx"))
    df = permanovas(taxprofiles[:, idx], [:stool_age, Symbol.(eeg_features)...])
    df.cohort .= label
    df.profile .= "taxa"
    append!(pnovs, df)
end

for (idx, label) in zip((v1_idx, v2_idx, v3_idx, v1v2_idx, v1v3_idx, v2v3_idx),
                        ("v1_idx", "v2_idx", "v3_idx", "v1v2_idx", "v1v3_idx", "v2v3_idx"))
    df = permanovas(unirefs[:, idx], [:stool_age, Symbol.(eeg_features)...])
    df.cohort .= label
    df.profile .= "unirefs"
    append!(pnovs, df)
end

using MultipleTesting
transform!(pnovs, "pvalue"=> (p-> adjust(p, BenjaminiHochberg()))=> "qvalue")
