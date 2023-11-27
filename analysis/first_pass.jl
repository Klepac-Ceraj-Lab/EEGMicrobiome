using FeatureSetEnrichments
using VKCComputing
using EEGMicrobiome

eeg = load_eeg()
mbo = load_microbiome(eeg.subject)


