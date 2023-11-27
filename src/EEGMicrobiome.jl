module EEGMicrobiome

export load_eeg,
       load_microbiome,
       load_taxonomic_profiles

using CSV
using DataFrames
using Chain
using VKCComputing

const _subject_excludes = Set([
    "191-29736014", # premature
    "191-34303377", # premature
])

function load_eeg(datadir = "./data")
    df = DataFrame()
    for tp in ("3m", "6m", "12m")
        path = joinpath(datadir, "eeg_vep_khulaSA_$tp.csv")
        tpdf = @chain path begin
            CSV.read(DataFrame, stringtype=String)
            subset!("infant_eye_anomaly_1yes_0no" => ByRow(x-> x != 1), 
                    "subject_id"=> ByRow(x-> x ∉ _subject_excludes)
            )
            select!("subject_id"               => "subject",
                    "age_$(tp)_vep" => "age", 
                    "Number_Segs_Post-Seg_Rej" => "n_segments",
                    r"^peak")
            rename(col-> replace(col, "_$tp"=> ""), _)
        end
        tpdf.timepoint .= tp
        append!(df, tpdf)
    end
    return df
end

function load_microbiome(subjects)
    base = LocalBase()
    df = DataFrame(subject = subjects, subject_rec = [get(base["Subjects"], "khula-$subject", missing) for subject in subjects])
    df.biospecimen_ids = [ismissing(rec) ? [missing] : get(rec, :Biospecimens, [missing]) for rec in df.subject_rec]
    df = flatten(df, :biospecimen_ids)
    df.biospecimen = [ismissing(rec) ? missing : base[rec][:uid] for rec in df.biospecimen_ids]
    df.seqprep_ids = [ismissing(rec) ? [missing] : get(base[rec], :seqprep, [missing]) for rec in df.biospecimen_ids]
    df = flatten(df, :seqprep_ids)
    subset!(df, :seqprep_ids=> ByRow(s-> !ismissing(s) && base[s][:keep] == 1)) 
    df.seqprep = [base[rec][:uid] for rec in df.seqprep_ids]
    df.age = Union{Missing, Float64}[get(base[rec], :subject_age, missing) for rec in df.biospecimen_ids]
    return df
    
end

function load_taxonomic_profiles(seqs)
end

end