const _subject_excludes = Set([
    "191-29736014", # premature
    "191-34303377", # premature
])

function load_cohorts(tab = "./data/allmeta.csv")
    df = CSV.read(tab, DataFrame; stringtype=String)
    subset!(df, "subject_id"=> ByRow(s-> s ∉ _subject_excludes))
    return df
end

function get_cohort(metadf, cohort)
    m = match(r"^(v\d)(v\d)?$", cohort)
    isnothing(m) && throw(ArgumentError("$cohort is not a valid cohort, use \"v1\", \"v2v3\" etc"))
    if isnothing(m[2])
        df = subset(metadf, "cohort_$cohort"=> identity)
        transform!(df, "age_vep_weeks" => ByRow(aw-> (aw / 52) * 12) => "vep_age")
        transform!(df, AsTable(["stool_age", "vep_age"]) => ByRow(nt-> nt.stool_age - nt.vep_age) => "age_diff")
        return df
    else
        subdf = subset(metadf, AsTable(Regex("cohort_$(cohort)_")) => ByRow(any))
        grp = groupby(subdf, "subject_id")
        @assert all(==(2), combine(grp, "visit"=> length => "n").n)
        @assert all(sub -> !ismissing(sub.stool_age[1]) && !ismissing(sub.age_vep_weeks[2]), grp)
        @assert all(sub -> issorted(sub.visit), grp)

        subdf = select(combine(grp,
            AsTable(["subject_id", "stool_age", "biospecimen", "seqprep", "S_well", "filename", "taxprofile", "genefamilies"]) => 
                (sub-> NamedTuple(Symbol(c) => first(sub[Symbol(c)]) for c in keys(sub))) => AsTable,
            AsTable(r"peak|n_trial")=> (sub-> NamedTuple(Symbol(c) => last(sub[Symbol(c)]) for c in keys(sub))) => AsTable,
            "age_vep_weeks" => (aw-> (last(aw) / 52) * 12) => "vep_age",
            "visit" => (v-> cohort) => "visit"
            ), "subject_id", "visit", "stool_age", "vep_age", "biospecimen", "n_trials", "seqprep", r"peak_", "S_well", "filename", "taxprofile", "genefamilies"
        )
        subdf.age_diff = subdf.stool_age .- subdf.vep_age
        return subdf
    end
end
