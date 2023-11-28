using FeatureSetEnrichments
using VKCComputing
using EEGMicrobiome

# load data
eeg = load_eeg()
mbo = load_microbiome(eeg.subject)

load_taxonomic_profiles!(mbo)
load_functional_profiles!(mbo)
unique!(mbo, :seqprep)
sort!(mbo, "age")
subset!(mbo, "age"=> ByRow(!ismissing))
# 3M concurrent EEG

eeg3m = subset(eeg, "timepoint"=> ByRow(==("3m")))
leftjoin!(eeg3m, select(unique(subset(DataFrame(get(taxa)), "age"=> ByRow(a-> a < 9)), :subject), "subject", "age"=>"stool_age", "sample"); on=:subject)
subset!(eeg3m, "sample"=> ByRow(!ismissing))

leftjoin!(eeg3m, unique(select(subset(mbo, "age"=>ByRow(<(6))), Not("age")), "subject"); on=:subject)


