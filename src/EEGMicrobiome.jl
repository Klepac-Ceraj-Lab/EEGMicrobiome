module EEGMicrobiome

export load_eeg,
       load_microbiome,
       load_taxonomic_profiles!,
       load_functional_profiles!

using CSV
using DataFrames
using Chain
using VKCComputing
using Preferences
using BiobakeryUtils: metaphlan_profiles
using ThreadsX
using HypothesisTests: MannWhitneyUTest, pvalue
using MultipleTesting: adjust, BenjaminiHochberg
using GLM


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
            transform!("age"=> ByRow((a-> a / 52 * 12))=> "age") # convert to months
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
    df.visit = [ismissing(first(get(base[rec], :visit, [missing]))) ? missing : base[first(base[rec][:visit])][:uid] for rec in df.biospecimen_ids]
    return select(df, "subject", "biospecimen", "seqprep", "age", "visit")
    
end

function load_taxonomic_profiles!(mbiome_tab)
    seqs = mbiome_tab.seqprep

    mgxdir = load_preference(VKCComputing, "mgx_analysis_dir")
    files = filter(readdir(joinpath(mgxdir, "metaphlan"); join=true)) do file
        m = match(r"^SEQ\d+", basename(file))
        !isnothing(m) && endswith(basename(file), "_profile.tsv") && m.match ∈ seqs
    end
    taxa = mapreduce(vcat, files) do f
        seqprep = replace(basename(f), r"_S\d+_profile\.tsv"=> "")
        df = CSV.read(f, DataFrame; skipto = 5, header=["feature", "taxid", "abundance", "addtl_spec"])
        df.seqprep .= seqprep
        subset!(df, "feature"=> ByRow(f-> contains(f, "|s__")))
        transform!(df, "feature"=> ByRow(f-> last(split(f, '|')))=>"feature")
        select!(df, "seqprep", "feature", "abundance")
    end
    taxa_wide = unstack(taxa, "feature", "abundance")
    foreach(names(taxa_wide)) do n
        taxa_wide[!, n] = coalesce.(taxa_wide[!, n], 0.)
    end
    leftjoin!(mbiome_tab, taxa_wide; on=:seqprep)

end

function load_functional_profiles!(mbiome_tab)
    seqs = mbiome_tab.seqprep

    mgxdir = load_preference(VKCComputing, "mgx_analysis_dir")
    files = filter(readdir(joinpath(mgxdir, "humann", "main"); join=true)) do file
        m = match(r"^SEQ\d+", basename(file))
        !isnothing(m) && endswith(basename(file), "_genefamilies.tsv") && m.match ∈ seqs
    end
    gfs = mapreduce(vcat, files) do f
        seqprep = replace(basename(f), r"_S\d+_genefamilies\.tsv"=> "")
        df = CSV.read(f, DataFrame; skipto = 2, header=["feature", "abundance"])
        df.seqprep .= seqprep
        subset!(df, "feature"=> ByRow(f-> !contains(f, '|')))
        select!(df, "seqprep", "feature", "abundance")
    end
    gfs_wide = unstack(gfs, "feature", "abundance")
    foreach(names(gfs_wide)) do n
        gfs_wide[!, n] = coalesce.(gfs_wide[!, n], 0.)
    end
    leftjoin!(mbiome_tab, gfs_wide; on=:seqprep)
end


function runlms(indf, outfile, respcol, featurecols;
                age_col = "age",
                additional_cols = Symbol[],
                formula = term(:func) ~ term(respcol) + term(age_col) + term(:n_segments),
                prevalence_filter = 0.05
        )
        @debug "Respcol: $respcol"
    lmresults = DataFrame(ThreadsX.map(featurecols) do feature
        # @debug "Feature: $feature"

        # ab = collect(indf[!, feature] .+ (minimum(indf[over0, feature])) / 2) # add half-minimum non-zerovalue
        df = select(indf, respcol, age_col, "n_segments", additional_cols...)
        over0 = indf[!, feature] .> 0
        default_ret = (; feature, Name = respcol, coef = NaN, std_err = NaN, z = NaN, pvalue = NaN, lower_95 = NaN, upper_95 = NaN, qvalue = NaN)
        
        (sum(skipmissing(over0)) / length(over0)) > prevalence_filter || return default_ret

        df.func = over0
        # @debug "DataFrame: $df"

        try
            mod = glm(formula, df, Binomial(), LogitLink(); dropcollinear=false)
            ct = DataFrame(coeftable(mod))
            ct.feature .= feature
            rename!(ct, "Pr(>|z|)"=>"pvalue", "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
            select!(ct, Cols(:feature, :Name, :))
            return NamedTuple(only(filter(row-> row.Name == respcol, eachrow(ct))))    
        catch e
            @warn "hit $e for $feature"
            return default_ret
        end
    end)

    subset!(lmresults, "pvalue"=> ByRow(!isnan))
    DataFrames.transform!(lmresults, :pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(lmresults, :qvalue)

    CSV.write(outfile, lmresults)
    lmresults
end


end