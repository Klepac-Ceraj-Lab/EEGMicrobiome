module EEGMicrobiome

export load_eeg

using CSV
using DataFrames
using Chain

const _subject_excludes = Set([
    "191-29736014", # premature
    "191-34303377", # premature
])

function load_eeg(datadir = "./data")
    df = DataFrame()
    for tp in ("3m", "6m", "12m")
        path = joinpath(datadir, "eeg_vep_khulaSA_$tp.csv")
        tpdf = @chain path begin
            CSV.read(DataFrame)
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

end