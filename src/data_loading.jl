const _subject_excludes = Set([
    "191-29736014", # premature
    "191-34303377", # premature
])

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

function load_cohort(cohort)
    tab = CSV.read("data/allmeta.csv", DataFrame)
    subset!(tab, "visit"=> ByRow(==(cohort)))
    tab.peak_latency_P1_corrected = tab.peak_latency_P1 .- tab.peak_latency_N1
    tab.peak_latency_N2_corrected = tab.peak_latency_N2 .- tab.peak_latency_P1
    tab.peak_amp_P1_corrected = tab.peak_amp_P1 .- tab.peak_amp_N1
    tab.peak_amp_N2_corrected = tab.peak_amp_N2 .- tab.peak_amp_P1
 

    # length(files) == size(tab, 1) || throw(ArgumentError(
    #     "Mismatch between eeg and microbiome data, missing: $(
    #         setdiff(tab.seqprep, map(f-> replace(basename(f), r"_S\d+.+"=>""), files)))"
    #    ))
    comm = metaphlan_profiles(files)
    comm = CommunityProfile(abundances(comm), features(comm), MicrobiomeSample.(replace.(samplenames(comm), r"_S\d+_profile"=>""))) 

    set!(comm, select(tab, "seqprep"=>"sample", Cols(:)))
    return comm
end
