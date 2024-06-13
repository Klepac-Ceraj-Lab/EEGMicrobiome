using CSV
using DataFrames
using Chain
using VKCComputing

vep = CSV.read("data/allvep.csv", DataFrame; stringtype=String)
demo = filter(row-> !all(ismissing, row), CSV.read("data/khula_metadata.csv", DataFrame; stringtype=String))

#-


base = LocalBase(; update=true)

airtab_subject = mapreduce((df1, df2)-> vcat(df1, df2; cols=:union), vep.subject_id) do sub
    subuid = "khula-$sub"
    newrec = (; subject_id = sub, subject_uid = subuid)
    rec = get(base["Subjects"], subuid, missing)

    return ismissing(rec) ? DataFrame([newrec]) :
			    DataFrame([(; newrec..., rec.fields...)])
end

dropmissing!(airtab_subject)
subset!(airtab_subject, "keep" => ByRow(==(1)))
airtab_subject = flatten(airtab_subject, "Biospecimens")
airtab_subject.project = only.(airtab_subject.project)

airtab_biospecimens = mapreduce((df1, df2)-> vcat(df1, df2; cols=:union), eachrow(airtab_subject)) do row    
    biosp = row.Biospecimens
    newrec = (; biospecimen = biosp, subject_id = row.subject_id)
    rec = get(base["Biospecimens"], biosp, missing)

    return ismissing(rec) ? DataFrame([newrec]) :
			    DataFrame([(; newrec..., rec.fields...)])
end

subset!(airtab_biospecimens, "keep" => ByRow(==(1)), "seqprep"=> ByRow(!ismissing))
airtab_biospecimens = flatten(airtab_biospecimens, ["seqprep", "sequencing_batch (from seqprep)"])
dropmissing!(airtab_biospecimens, "visit")
airtab_biospecimens.visit = [base[v].fields.uid for v in only.(airtab_biospecimens.visit)]

airtab_seqpreps = mapreduce((df1, df2)-> vcat(df1, df2; cols=:union), eachrow(airtab_biospecimens)) do row    
    seqprep = row.seqprep
    biospecimen = row.uid
    subject_id = row.subject_id
    visit = row.visit
    rec = get(base["SequencingPrep"], seqprep, missing)

    if ismissing(rec)
	return DataFrame([(; subject_id, visit, biospecimen, seqprep)])
    else
    	fields = rec.fields
	seqprep = fields.uid
	return DataFrame([(; subject_id, visit, biospecimen, seqprep,
			   NamedTuple(k=> get(fields, k, missing) for k in [:S_well, :filename, :kneaddata, :metaphlan, :humann, :keep])...)])
    end
end

subset!(airtab_seqpreps, "keep"=> ByRow(==(1)))
dropmissing!(airtab_seqpreps)
unique!(airtab_seqpreps)

@assert airtab_seqpreps
