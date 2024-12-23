# start from Figures.jl - get to `mdata`

nofile = filter(!isfile, String.(skipmissing(mdata.taxprofile)))
df = DataFrame(profile=basename.(nofile));

df.seq = replace.(basename.(nofile), "_profile.tsv"=> "")
kneaddata = DataFrame(
    path = filter(f-> replace(basename(f), r"(SEQ\d+_S\d+).+" => s"\1") ∈ df.seq &&
                      contains(basename(f), r"kneaddata_(unmatched|paired)"),
                  readdir("/grace/sequencing/processed/mgx/kneaddata/"; join=true))
)

transform!(kneaddata, "path"=> ByRow(path-> begin
    file = basename(path)
    dir = dirname(path)
    seq = replace(file, r"(SEQ\d+_S\d+).+"=> s"\1")
    (; file, dir, seq)
end)=> ["file", "dir", "seq"])

leftjoin!(kneaddata, df; on="seq")
kneaddata.profile = [p for p in kneaddata.profile]

grouped = groupby(kneaddata, "seq")

sem = Base.Semaphore(2)

@sync for seq in grouped
    Threads.@spawn begin
        Base.acquire(sem) do
            seqname = first(seq.seq)
            @info first(seq.seq)
            @info "    concatenating fastq files"
            knead_cat = "/murray/kevin/mp3_fuckup/kneaddata_cat/$(seqname)_cat.fastq.gz"
            isfile(knead_cat) || run(pipeline(`cat $(seq.path)`; stdout=knead_cat))

            @info "    running metaphlan"
            if !isfile("/murray/kevin/mp3_fuckup/metaphlan/$(seqname)_profile.tsv")
                cmd = Cmd([
                    "metaphlan", knead_cat,
                    "/murray/kevin/mp3_fuckup/metaphlan/$(seqname)_profile.tsv",
                    "--input_type", "fastq",
                    "--bowtie2out", "/murray/kevin/mp3_fuckup/metaphlan/$(seqname)_bowtie2.tsv",
                    "--samout", "/murray/kevin/mp3_fuckup/metaphlan/$(seqname).sam",
                    "--nproc", "16",
                    "--bowtie2db", "/murray/databases/metaphlan/"
                ])
                run(cmd)
            end
        end
    end
end
