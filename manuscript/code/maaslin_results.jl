using CSV
using DataFrames
using Chain

all_results = mapreduce(vcat, walkdir("data/outputs/maaslin/")) do (root, dirs, files)
    contains(basename(root), "model_") || return DataFrame()
    model = replace(basename(root), "model_"=>"")
    res = CSV.read(joinpath(root, "all_results.tsv"), DataFrame; delim='\t')
    res.model .= model
    return res
end

sort!(all_results, "qval_joint")

first(select(all_results, :feature, :model, :metadata, :qval_joint), 2)

@chain all_results begin
    subset("metadata"=> ByRow(m->!contains(m, "age")))
    select("model", "feature", "metadata", Not(["name","value"]))
    CSV.write("data/outputs/Table_maaslin_species.csv",_)
end
