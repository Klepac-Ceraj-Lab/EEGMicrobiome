using CSV
using DataFrames
using Chain
using VKCComputing
using Preferences

vep = CSV.read("data/allvep.csv", DataFrame; stringtype=String)
transform!(vep, "visit_type" => ByRow(v -> "$(v)mo") => "visit")
demo = filter(row -> !all(ismissing, row), CSV.read("data/khula_metadata.csv", DataFrame; stringtype=String))

#-


base = LocalBase(; update=true)

airtab_subject = mapreduce((df1, df2) -> vcat(df1, df2; cols=:union), vep.subject_id) do sub
    subuid = "khula-$sub"
    newrec = (; subject_id=sub, subject_uid=subuid)
    rec = get(base["Subjects"], subuid, missing)

    return ismissing(rec) ? DataFrame([newrec]) :
           DataFrame([(; newrec..., rec.fields...)])
end

dropmissing!(airtab_subject)
subset!(airtab_subject, "keep" => ByRow(==(1)))
airtab_subject = flatten(airtab_subject, "Biospecimens")
airtab_subject.project = only.(airtab_subject.project)

airtab_biospecimens = mapreduce((df1, df2) -> vcat(df1, df2; cols=:union), eachrow(airtab_subject)) do row
    biosp = row.Biospecimens
    newrec = (; biospecimen=biosp, subject_id=row.subject_id)
    rec = get(base["Biospecimens"], biosp, missing)

    return ismissing(rec) ? DataFrame([newrec]) :
           DataFrame([(; newrec..., rec.fields...)])
end

subset!(airtab_biospecimens, "keep" => ByRow(==(1)), "seqprep" => ByRow(!ismissing))
airtab_biospecimens = flatten(airtab_biospecimens, ["seqprep", "sequencing_batch (from seqprep)"])
dropmissing!(airtab_biospecimens, "visit")
airtab_biospecimens.visit = [base[v].fields.uid for v in only.(airtab_biospecimens.visit)]

airtab_seqpreps = mapreduce((df1, df2) -> vcat(df1, df2; cols=:union), eachrow(airtab_biospecimens)) do row
    seqprep = row.seqprep
    biospecimen = row.uid
    subject_id = row.subject_id
    visit = row.visit
    stool_age = row.subject_age
    rec = get(base["SequencingPrep"], seqprep, missing)

    if ismissing(rec)
        return DataFrame([(; subject_id, visit, biospecimen, seqprep, stool_age)])
    else
        fields = rec.fields
        seqprep = fields.uid
        return DataFrame([(; subject_id, visit, biospecimen, seqprep, stool_age,
            NamedTuple(k => get(fields, k, missing) for k in [:S_well, :filename, :kneaddata, :metaphlan, :humann, :keep])...)])
    end
end

subset!(airtab_seqpreps, "keep" => ByRow(==(1)))
dropmissing!(airtab_seqpreps)
unique!(airtab_seqpreps)

@assert nrow(airtab_seqpreps) == count(airtab_seqpreps.kneaddata)
@assert nrow(airtab_seqpreps) == count(airtab_seqpreps.metaphlan)
@assert nrow(airtab_seqpreps) == count(airtab_seqpreps.humann)


analysis_dir = load_preference(VKCComputing, "mgx_analysis_dir", ".")

transform!(airtab_seqpreps, "filename" => ByRow(filename -> begin
    taxprofile = joinpath(analysis_dir, "metaphlan", "$(filename)_profile.tsv")
    genefamilies = joinpath(analysis_dir, "humann", "main", "$(filename)_genefamilies.tsv")

    (; taxprofile, genefamilies)
end) => ["taxprofile", "genefamilies"]
)

@assert all(isfile, airtab_seqpreps.taxprofile)
@assert all(isfile, airtab_seqpreps.genefamilies)

allmeta = outerjoin(vep, airtab_seqpreps; on=["subject_id", "visit"])
rename!(allmeta, "Number_Segs_Post.Seg_Rej" => "n_trials")
select!(allmeta, Cols(
    "subject_id",
    "visit",
    "age_vep_weeks",
    "stool_age",
    "biospecimen",
    "seqprep",
    "S_well",
    r"^peak_",
    "n_trials",
    "filename",
    "taxprofile",
    "genefamilies"
)
)


subset!(groupby(allmeta, "subject_id"), "stool_age" => (a -> !all(ismissing, a)))


function visitlt(v1, v2)
    v1 = parse(Int, replace(v1, "mo" => ""))
    v2 = parse(Int, replace(v2, "mo" => ""))
    v1 < v2
end

allmeta = @chain allmeta begin
    sort("visit"; lt=visitlt)
    groupby("subject_id")
    subset(AsTable(["age_vep_weeks", "stool_age"]) => (nt -> begin
        va = nt.age_vep_weeks
        sa = nt.stool_age

        vi = findall(!ismissing, va)
        si = findall(!ismissing, sa)

        any(s -> any(v -> v >= s, vi), si)
    end))
end

transform!(groupby(allmeta, "subject_id"), AsTable(["visit", "age_vep_weeks", "stool_age"]) => (
    nt -> begin
        sidx = findall(!ismissing, nt.stool_age)
        vidx = findall(!ismissing, nt.age_vep_weeks)

        cohort_v1 = map(i -> nt.visit[i] == "3mo" && i in sidx && i in vidx, eachindex(nt.visit))
        cohort_v2 = map(i -> nt.visit[i] == "6mo" && i in sidx && i in vidx, eachindex(nt.visit))
        cohort_v3 = map(i -> nt.visit[i] == "12mo" && i in sidx && i in vidx, eachindex(nt.visit))


        cohort_v1v2_stool = map(i -> nt.visit[i] == "3mo" && i in sidx && any(j -> nt.visit[j] == "6mo", vidx), eachindex(nt.visit))
        cohort_v1v3_stool = map(i -> nt.visit[i] == "3mo" && i in sidx && any(j -> nt.visit[j] == "12mo", vidx), eachindex(nt.visit))
        cohort_v2v3_stool = map(i -> nt.visit[i] == "6mo" && i in sidx && any(j -> nt.visit[j] == "12mo", vidx), eachindex(nt.visit))

        cohort_v1v2_vep = map(i -> nt.visit[i] == "6mo" && i in vidx && any(j -> nt.visit[j] == "3mo", sidx), eachindex(nt.visit))
        cohort_v1v3_vep = map(i -> nt.visit[i] == "12mo" && i in vidx && any(j -> nt.visit[j] == "3mo", sidx), eachindex(nt.visit))
        cohort_v2v3_vep = map(i -> nt.visit[i] == "12mo" && i in vidx && any(j -> nt.visit[j] == "6mo", sidx), eachindex(nt.visit))

        cohort_v1v1v2_stool = map(i -> nt.visit[i] == "3mo" && i in sidx &&
                                           any(j -> nt.visit[j] == "3mo", vidx) &&
                                           any(j -> nt.visit[j] == "6mo", vidx), eachindex(nt.visit))
        cohort_v1v1v3_stool = map(i -> nt.visit[i] == "3mo" && i in sidx &&
                                           any(j -> nt.visit[j] == "3mo", vidx) &&
                                           any(j -> nt.visit[j] == "12mo", vidx), eachindex(nt.visit))
        cohort_v2v2v3_stool = map(i -> nt.visit[i] == "6mo" && i in sidx &&
                                           any(j -> nt.visit[j] == "6mo", vidx) &&
                                           any(j -> nt.visit[j] == "12mo", vidx), eachindex(nt.visit))

        # cohorts that have prior stool AND prior vep
        cohort_v1v1v2_vep = map(i -> nt.visit[i] == "3mo" && i in vidx &&
                                         any(j -> nt.visit[j] == "3mo", sidx) &&
                                         any(j -> nt.visit[j] == "6mo", vidx), eachindex(nt.visit))
        cohort_v1v1v3_vep = map(i -> nt.visit[i] == "3mo" && i in vidx &&
                                         any(j -> nt.visit[j] == "3mo", sidx) &&
                                         any(j -> nt.visit[j] == "12mo", vidx), eachindex(nt.visit))
        cohort_v2v2v3_vep = map(i -> nt.visit[i] == "6mo" && i in vidx &&
                                         any(j -> nt.visit[j] == "6mo", sidx) &&
                                         any(j -> nt.visit[j] == "12mo", vidx), eachindex(nt.visit))

        # Cohorts with vep prior to stool (for control)
        cohort_v2v1_stool = map(i -> nt.visit[i] == "6mo" && i in sidx && any(j -> nt.visit[j] == "3mo", vidx), eachindex(nt.visit))
        cohort_v3v1_stool = map(i -> nt.visit[i] == "12mo" && i in sidx && any(j -> nt.visit[j] == "3mo", vidx), eachindex(nt.visit))
        cohort_v3v2_stool = map(i -> nt.visit[i] == "12mo" && i in sidx && any(j -> nt.visit[j] == "6mo", vidx), eachindex(nt.visit))

        cohort_v2v1_vep = map(i -> nt.visit[i] == "3mo" && i in vidx && any(j -> nt.visit[j] == "6mo", sidx), eachindex(nt.visit))
        cohort_v3v1_vep = map(i -> nt.visit[i] == "3mo" && i in vidx && any(j -> nt.visit[j] == "12mo", sidx), eachindex(nt.visit))
        cohort_v3v2_vep = map(i -> nt.visit[i] == "6mo" && i in vidx && any(j -> nt.visit[j] == "12mo", sidx), eachindex(nt.visit))

        return (; cohort_v1, cohort_v2, cohort_v3,
            cohort_v1v2_stool, cohort_v1v3_stool, cohort_v2v3_stool,
            cohort_v1v2_vep, cohort_v1v3_vep, cohort_v2v3_vep,
            cohort_v1v1v2_stool, cohort_v1v1v3_stool, cohort_v2v2v3_stool,
            cohort_v1v1v2_vep, cohort_v1v1v3_vep, cohort_v2v2v3_vep,
            cohort_v2v1_stool, cohort_v3v1_stool, cohort_v3v2_stool,
            cohort_v2v1_vep, cohort_v3v1_vep, cohort_v3v2_vep,
        )
    end
) => AsTable
)

subset!(allmeta, AsTable(r"cohort_") => ByRow(nt -> any(!ismissing, nt)))
transform!(allmeta,
    AsTable(["peak_latency_P1", "peak_latency_N1"]) => (tab -> tab.peak_latency_P1 .- tab.peak_latency_N1) => "peak_latency_P1_corrected",
    AsTable(["peak_amp_P1", "peak_amp_N1"]) => (tab -> tab.peak_amp_P1 .- tab.peak_amp_N1) => "peak_amp_P1_corrected",
    AsTable(["peak_latency_N2", "peak_latency_P1"]) => (tab -> tab.peak_latency_N2 .- tab.peak_latency_P1) => "peak_latency_N2_corrected",
    AsTable(["peak_amp_N2", "peak_amp_P1"]) => (tab -> tab.peak_amp_N2 .- tab.peak_amp_P1) => "peak_amp_N2_corrected",
    "visit" => ByRow(v -> replace(v, "3mo" => "v1", "6mo" => "v2", "12mo" => "v3")) => "visit"
)



CSV.write("data/allmeta.csv", allmeta)



