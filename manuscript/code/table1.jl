using VKCComputing
using EEGMicrobiome
using DataFrames
using Chain

# load data
eeg = load_eeg()
eeg.peak_latency_P1_corrected = eeg.peak_latency_P1 .- eeg.peak_latency_N1
eeg.peak_latency_N2_corrected = eeg.peak_latency_N2 .- eeg.peak_latency_P1
eeg.peak_amp_P1_corrected = eeg.peak_amp_P1 .- eeg.peak_amp_N1
eeg.peak_amp_N2_corrected = eeg.peak_amp_N2 .- eeg.peak_amp_P1
rename!(eeg, "age"=> "eeg_age")
mbo = load_microbiome(eeg.subject)
transform!(mbo, "visit"=>ByRow(v-> replace(v, "mo"=>"m"))=> "visit")

eegmbo = @chain eeg begin
    select("subject", "timepoint"=>"visit", Cols(:))
    unique!(["subject", "timepoint"])
    outerjoin(mbo; on = ["subject", "visit"])
    groupby("subject")
    subset(AsTable(["visit", "seqprep", "eeg_age"])=> (nt-> begin
	# skip subjects without at least 1 eeg and at least 1 seqprep
	(all(ismissing, nt.seqprep) || all(ismissing(nt.eeg_age))) && return false
	# keep subjects with at least 1 concurrent eeg/seqprep
	any(.!ismissing.(nt.seqprep) .& .!ismissing.(nt.eeg_age)) && return true
	# Keep any remaining subjects that have a seqprep *prior* to an EEG
	srt = sortperm([parse(Int, replace(v, "m"=>"")) for v in nt.visit])
	any(findfirst(!ismissing, nt.seqprep[srt]) .<= findfirst(!ismissing, nt.eeg_age[srt]))
    end))
end

wide_sub = select(leftjoin(
    select(unstack(eegmbo, "subject", "visit", "eeg_age"), "subject", "3m"=>"eeg_3m", "6m"=> "eeg_6m", "12m"=>"eeg_12m"),
    select(unstack(eegmbo, "subject", "visit", "age"), "subject", "3m"=>"seqprep_3m", "6m"=> "seqprep_6m", "12m"=>"seqprep_12m"),
    on="subject"), "subject", r"3m", r"6m", r"12m")

transform!(wide_sub, AsTable(r"3m")=> ByRow(nt-> nt.eeg_3m - nt.seqprep_3m) => "diff_3m",
		    AsTable(r"6m")=> ByRow(nt-> nt.eeg_6m - nt.seqprep_6m) => "diff_6m",
		   AsTable(r"12m")=> ByRow(nt-> nt.eeg_12m - nt.seqprep_12m) => "diff_12m")

CSV.write("data/outputs/eeg_microbiome_subjects.csv", wide_sub)
CSV.write("data/outputs/eeg_microbiome_timepoints.csv", sort(select(eegmbo, "subject", "visit", "eeg_age", "seqprep"), ["subject", "visit"]))

# Generate table subsets
concurrent_3m = select(subset(wide_sub, AsTable(r"3m") => ByRow((nt-> all(!ismissing, values(nt))))), "subject")
concurrent_6m = select(subset(wide_sub, AsTable(r"6m") => ByRow((nt-> all(!ismissing, values(nt))))), "subject")
concurrent_12m = select(subset(wide_sub, AsTable(r"12m") => ByRow((nt-> all(!ismissing, values(nt))))), "subject")

future_3m6m = select(subset(wide_sub, AsTable(["seqprep_3m", "eeg_6m"]) => ByRow((nt-> all(!ismissing, values(nt))))), "subject")
future_3m12m = select(subset(wide_sub, AsTable(["seqprep_3m", "eeg_12m"]) => ByRow((nt-> all(!ismissing, values(nt))))), "subject")
future_6m12m = select(subset(wide_sub, AsTable(["seqprep_6m", "eeg_12m"]) => ByRow((nt-> all(!ismissing, values(nt))))), "subject")
 
# join relevant data

dropmissing!(leftjoin!(concurrent_3m, select(subset(eegmbo, "visit"=> ByRow(==("3m"))), "subject", "seqprep", "age", "n_segments", r"peak_"); on="subject"))
dropmissing!(leftjoin!(concurrent_6m, select(subset(eegmbo, "visit"=> ByRow(==("6m"))), "subject", "seqprep", "age", "n_segments", r"peak_"); on="subject"))
dropmissing!(leftjoin!(concurrent_12m, select(subset(eegmbo, "visit"=> ByRow(==("12m"))), "subject", "seqprep", "age", "n_segments", r"peak_"); on="subject"))


@chain future_3m6m begin
	leftjoin!(select(subset(eegmbo, "visit"=>ByRow(==("3m"))), "subject", "seqprep", "age"); on="subject")
	leftjoin!(select(subset(eegmbo, "visit"=>ByRow(==("6m"))), "subject", "eeg_age", "n_segments", r"peak_"); on="subject")
	dropmissing!
end
@chain future_3m12m begin
	leftjoin!(select(subset(eegmbo, "visit"=>ByRow(==("3m"))), "subject", "seqprep", "age"); on="subject")
	leftjoin!(select(subset(eegmbo, "visit"=>ByRow(==("12m"))), "subject", "eeg_age", "n_segments", r"peak_"); on="subject")
	dropmissing!
end
@chain future_6m12m begin
	leftjoin!(select(subset(eegmbo, "visit"=>ByRow(==("6m"))), "subject", "seqprep", "age"); on="subject")
	leftjoin!(select(subset(eegmbo, "visit"=>ByRow(==("12m"))), "subject", "eeg_age", "n_segments", r"peak_"); on="subject")
	dropmissing!
end

cohortsdir = "data/outputs/cohort_tables"
isdir(cohortsdir) || mkdir(cohortsdir)

CSV.write(joinpath(cohortsdir, "concurrent_3m.csv"), concurrent_3m)
CSV.write(joinpath(cohortsdir, "concurrent_6m.csv"), concurrent_6m)
CSV.write(joinpath(cohortsdir, "concurrent_12m.csv"), concurrent_12m)
CSV.write(joinpath(cohortsdir, "future_3m6m.csv"), future_3m6m)
CSV.write(joinpath(cohortsdir, "future_3m12m.csv"), future_3m12m)
CSV.write(joinpath(cohortsdir, "future_6m12m.csv"), future_6m12m)

##

fig = Figure(; size=(800, 1200))
ax1 = Axis(fig[1, 1]; title = "3m", xlabel="eeg age - microbiome age")
ax2 = Axis(fig[2, 1]; title = "6m", xlabel="eeg age - microbiome age")
ax3 = Axis(fig[3, 1]; title = "12m", xlabel="eeg age - microbiome age")

hist!(ax1, collect(skipmissing(wide_sub.diff_3m)))
hist!(ax2, collect(skipmissing(wide_sub.diff_6m)))
hist!(ax3, collect(skipmissing(wide_sub.diff_12m)))

save("data/figures/age_diff_hist.png", current_figure())
