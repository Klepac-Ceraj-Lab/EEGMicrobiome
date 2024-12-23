const _subject_excludes = Set([
    "191-29736014", # premature
    "191-34303377", # premature
])

function load_cohorts(tab="./data/allmeta.csv")
    df = CSV.read(tab, DataFrame; stringtype=String)
    subset!(df, "subject_id" => ByRow(s -> s ∉ _subject_excludes))
    transform!(df, "age_vep_weeks" => ByRow(a -> a / 52 * 12) => "eeg_age")
    return df
end

function get_cohort(metadf, cohort)
    m = match(r"^(v\d)(v\d)?$", cohort)
    isnothing(m) && throw(ArgumentError("$cohort is not a valid cohort, use \"v1\", \"v2v3\" etc"))
    if isnothing(m[2])
        df = subset(metadf, "cohort_$cohort" => identity)
        transform!(df, "age_vep_weeks" => ByRow(aw -> (aw / 52) * 12) => "vep_age")
        transform!(df, AsTable(["stool_age", "vep_age"]) => ByRow(nt -> nt.stool_age - nt.vep_age) => "age_diff")
        return df
    else
        subdf = subset(metadf, AsTable(Regex("cohort_$(cohort)")) => ByRow(any))
        sort!(subdf, "subject_id")
        subdf_stool = subset(subdf, "cohort_$(cohort)_stool" => identity)
        subdf_vep = subset(subdf, "cohort_$(cohort)_vep" => identity)
        transform!(subdf_vep, "age_vep_weeks" => ByRow(aw -> (aw / 52) * 12) => "vep_age")

        @assert all(subdf_stool.subject_id .== subdf_vep.subject_id)
        @assert all(!ismissing, subdf_stool.stool_age)
        @assert all(!ismissing, subdf_vep.vep_age)

        subdf = leftjoin(
            select(subdf_stool, "subject_id", "stool_age", "biospecimen", "seqprep", "S_well", "filename", "taxprofile", "genefamilies"),
            select(subdf_vep, Not(["stool_age", "biospecimen", "seqprep", "S_well", "filename", "taxprofile", "genefamilies"]));
            on="subject_id"
        )
        @assert all(!ismissing, [subdf.stool_age; subdf.vep_age])
        subdf.visit .= cohort
        select!(subdf, "subject_id", "visit", "stool_age", "vep_age", "biospecimen", "n_trials", "seqprep", r"peak_", "S_well", "filename", "taxprofile", "genefamilies")
        subdf.age_diff = subdf.stool_age .- subdf.vep_age
        subdf.seqprep = collect(skipmissing(subdf.seqprep))
        return subdf
    end
end
