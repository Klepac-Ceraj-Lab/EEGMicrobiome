#set page(paper: "us-letter", numbering: "1 of 1")
#set text(8pt, font: "Liberation Sans", lang: "en")

//double spacing
#set text(top-edge: 0.7em, bottom-edge: -0.3em)
#set par(leading: 1em)

#set align(center)

//this is a test

// Only in local using latest development version of typst
#set par.line(numbering: "1")
#let no-paren-cite(label) = {
  show regex("\(|\)"): none
  cite(label, form: "normal") 
}

#show figure.caption: it => context [
  #set align(left)
  #set text(8pt)
  *#it.supplement~#it.counter.display()#it.separator*#it.body
]

// 89 / 200 char
= Co-development of gut microbial metabolism and visual neural circuitry over human infancy

Kevin S. Bonham#super([1,2])\*,
Emma T. Margolis#super([3])\*,
Guilherme Fahur Bottino#super([1]),
Ana Sobrino#super([3]),
Fadheela Patel#super([4]),
Shelley McCann#super([1]),
Michal R. Zieff#super([5]),
Marlie Miles#super([5]),
Donna Herr#super([5]),
Lauren Davel#super([5]),
Cara Bosco#super([3]),
Khula South Africa Data Collection Team#super([5]),
Curtis Huttenhower#super([6])
Nicolò Pini#super([7,8]),
Daniel C. Alexander#super([9]),
Derek K. Jones#super([10]),
Steve C. R. Williams#super([11]),
Dima Amso#super([12]),
Melissa Gladstone#super([13]),
William P. Fifer#super([7,8]),
Kirsten A. Donald#super([5,14]),
Laurel J. Gabard-Durnam#super([3])#sym.dagger,
Vanja Klepac-Ceraj#super([1])#sym.dagger,

#set align(left)

\*These authors contributed equally#linebreak()
#sym.dagger co-corresponding authors#linebreak()

#super([1])Department of Biological Sciences, Wellesley College, Wellesley, MA, USA#linebreak()
#super([1])Department of Medicine, Tufts Medical Center, Boston, MA, USA#linebreak()
#super([3])Department of Phychology, Northeastern University, Boston, MA, USA#linebreak()
#super([4])Division of Medical Microbiology, University of Cape Town, Cape Town, Western Cape, ZAF#linebreak()
#super([5])Department of Paediatrics and Child Health, University of Cape Town, Cape Town, Western Cape, ZAF#linebreak()
#super([6])Department of Biostatistics, Harvard T.H. Chan School of Public Health, Boston, MA, USA#linebreak()
#super([7])Department of Psychiatry, Columbia University, Irving Medical Center, New York City, NY, USA#linebreak()
#super([8])Division of Developmental Neuroscience, New York State Psychiatric Institute, New York City, NY, USA#linebreak()
#super([9])Centre for Medical Image Computing, Department of Computer Science,
    University College London, London, United Kingdom#linebreak()
#super([10])Cardiff University Brain Research Imaging Centre, Cardiff University, Cardiff, Wales, United Kingdom#linebreak()
#super([11])Department of Neuroimaging, King’s College London, United Kingdom#linebreak()
#super([12])Department of Psychology, Columbia University, New York City, NY, USA#linebreak()
#super([13])Department of Women and Children’s Health, Institute of Life Course and Medical Science,
    Alder Hey Children’s NHS Foundation Trust, University of Liverpool, Liverpool, United Kingdom#linebreak()
#super([14])Neuroscience Institute, University of Cape Town, Cape Town, Western Cape, ZAF

// 50 / 55 char
*Running title*: Infant gut microbiome and visual neurodevelopment

*Keywords*: visual cortex development, metagenome, infant gut microbiome, visual-evoked potential, neuroplasticity

=== ORCID numbers:
#table(columns: (1fr, 1fr, 1fr, 1fr),
    [Kevin S. Bonham], [0000-0003-3200-7533],
    [Emma T. Margolis], [0000-0002-2036-8078],
    [Guilherme Fahur Bottino], [0000-0003-1953-1576],
    [Ana Sobrino], [],
    [Fadheela Patel], [0000-0001-5177-7416],
    [Shelley McCann], [0000-0002-9753-7968],
    [Michal Zieff], [0000-0001-9352-9947],
    [Marlie Miles], [],
    [Donna Herr], [],
    [Lauren Davel], [],
    [Cara Bosco], [],
    [Curtis Huttenhower], [0000-0002-1110-0096],
    [Nicolo Pini], [0000-0002-0839-6033],
    [Daniel C. Alexander], [0000-0003-2439-350X],
    [Derek K. Jones], [0000-0003-4409-8049],
    [Steve C. R. Williams], [0000-0003-4299-1941],
    [Dima Amso], [0000-0001-6798-4698],
    [Melissa Gladstone], [0000-0002-2579-9301],
    [William P. Fifer], [0000-0002-6936-9303],
    [Kirsten A. Donald], [0000-0002-0276-9660],
    [Laurel J. Gabard-Durnam], [0000-0002-4564-8068],
    [Vanja Klepac-Ceraj], [0000-0001-5387-5706],
)

#pagebreak()

//
// Leave comments like this
#set par(first-line-indent: 2em, spacing: 0.65em)
#set text(10pt)


== Abstract

Infancy is a time of elevated neuroplasticity supporting rapid brain and sensory development.
The gut microbiome,
also undergoing extensive developmental changes in early life,
may influence brain development through metabolism of neuroactive compounds.
Here, we leverage longitudinal data from 194 infants across the first 18 months of life
to show that microbial genes encoding enzymes that metabolize molecules
playing a key role in modulating early neuroplasticity
are associated with visual cortical neurodevelopment,
measured by the Visual-Evoked Potential (VEP).
Neuroactive compounds included neurotransmitters GABA and glutamate,
the amino acid tryptophan,
and short-chain fatty acids involved in myelination,
including acetate and butyrate.
Microbial gene sets around 4 months of age
were strongly associated with the VEP from around 9 to 14 months of age
and showed more associations than concurrently measured gene sets,
suggesting microbial metabolism in early life may affect subsequent neural plasticity and development.

== Introduction

The gut microbiome in early life has potential long-term implications for brain and body health.
One important way this influence can occur is through interactions with the central nervous system as a “microbial-gut-brain axis”
@rheePrinciplesClinicalImplications2009
@collinsInterplayIntestinalMicrobiota2012
@brettMicrobiotaGutBrain2019
The metabolic potential of the microorganisms that inhabit the gut vastly exceeds that of human cells alone,
with microbial genes outnumbering host genes by a hundredfold
@gilbertCurrentUnderstandingHuman2018.
In particular,
gut microbes have the ability to metabolize and synthesize many neuroactive compounds
@valles-colomerNeuroactivePotentialHuman2019.
However, the physiological relevance of this in humans has been difficult to quantify,
particularly during initial neurological development in early life.

Extensive work in preclinical models suggests that these neuroactive compounds
can influence the brain through both direct and indirect pathways.
For example, major neurotransmitters (e.g., glutamate, γ-aminobutyric acid (GABA), serotonin, and dopamine)
are readily synthesized and degraded by intestinal microbes and can enter circulation
and pass the blood-brain barrier to influence central nervous function
@janikMagneticResonanceSpectroscopy2016
@ahmedMicrobiotaderivedMetabolitesDrivers2022
@matsumotoCerebralLowmolecularMetabolites2013
@kawaseGutMicrobiotaMice2017.
Glutamatergic/GABA-ergic signaling is critical for balancing the brain’s excitatory
and inhibitory neurotransmission levels,
and alterations in the bi-directional glutamatergic/GABA-ergic signaling
between the gut microbiome and brain are implicated in several physical and mental health conditions
@filpaRoleGlutamatergicNeurotransmission2016
@bajGlutamatergicSignalingMicrobiotaGutBrain2019.
Similarly, the gut and the microbiome are critical to the regulation of metabolism
for the neurotransmitters serotonin and dopamine,
particularly through the metabolism of dietary tryptophan
@agusGutMicrobiotaRegulation2018.
Moreover, short-chain fatty acids (SCFAs) produced by the gut microbiome
may impact the brain directly by modulating neurotrophic factors,
glial and microglial maturation and myelination, and neuroinflammation
@dalileRoleShortchainFatty2019
@ernyMicrobiotaderivedAcetateEnables2021.
Other indirect pathways for gut microbial influence on the brain include vagus nerve stimulation,
neuroendocrine modulation and immune system regulation
@rheePrinciplesClinicalImplications2009.

Rapidly growing literature connects the metabolic potential of the gut microbiome
and brain function in humans
(reviewed in 
#no-paren-cite(<ahmedMicrobiotaderivedMetabolitesDrivers2022>)
#no-paren-cite(<parkerGutMicrobesMetabolites2019>)
#no-paren-cite(<aburtoGastrointestinalBrainBarriers2024>)),
but the overwhelming majority of this research is conducted in cohorts of adult participants.
Importantly, both the gut microbiome and the brain undergo dramatic and rapid development over the first postnatal years
@bonhamGutresidentMicroorganismsTheir2023
@lippeElectrophysiologicalMarkersVisuocortical2007
@deanModelingHealthyMale2014.
However, very little is currently known about how gut-brain influences emerge
or change during this critical window
@ahrensInfantMicrobesMetabolites2024a
@meyerAssociationGutMicrobiota2022
@callaghanMindGutAssociations2020.
Interrogating this early co-development in humans is, therefore,
key to both understanding adaptive gut-brain function and behavior and to informing strategies to support it.
Specifically, the visual cortex has been shown to be sensitive to gut microbiome modulations in adults
@canipeDiversityGutmicrobiomeRelated2021a
and in rodents
@luporiGutMicrobiotaEnvironmentally2022a,
yet the visual cortex undergoes its most rapid period of plasticity and maturation over infancy at the same time the microbiome changes most significantly
@deenOrganizationHighlevelVisual2017
@ellisRetinotopicOrganizationVisual2021
@kiorpesVisualDevelopmentPrimates2015.
Visual cortical functional maturation can be robustly indexed via electroencephalography (EEG)
with the Visual-Evoked Potential (VEP) response to visual stimuli from birth.
The VEP is especially useful for indexing visual neurodevelopment as its morphology
includes amplitude deflections as well as latencies to those deflections
reflecting maturation of function and structure, respectively.

Here, we investigated the longitudinal co-development of microbial metabolic potential quantified
via genes encoding enzymes that metabolize neuroactive compounds and visual neurodevelopment
as indexed by the VEP in a longitudinal community sample of 194 infants from
Gugulethu in Cape Town, South Africa,
recruited as part of the prospective longitudinal “Khula” Study
@zieffCharacterizingDevelopingExecutive2024.
Stool samples and EEG were each collected at up to 3 visits in the first 18 months of life.
Shotgun metagenomic sequencing was used to obtain microbial gene sequences from infant stool samples.
To index visual cortical functional development,
latencies and peak amplitudes were extracted from each component of the VEP
(i.e.,
first negative-going deflection, N1;
first positive-going deflection, P1; and
second negative-going deflection, N2),
producing six VEP features of interest.
We evaluated the concurrent association between microbial genes and VEP amplitudes and latencies,
and we tested prospective influences of microbial genes from early visits on VEP changes at later visits.
In this way, we were able to reveal the temporal dynamics of gut-brain co-development
within individuals during this most critical window of plasticity in both systems.

== Methods and Materials

=== Cohort

==== Participants and Study Design

Infants were recruited from local community clinics in Gugulethu,
an informal settlement in Cape Town, South Africa,
as part of a prospective longitudinal study
(most enrollments happened prenatally with 16% of infants enrolled shortly after birth;
#no-paren-cite(<zieffCharacterizingDevelopingExecutive2024>)).
The first language for the majority of residents in this area is Xhosa.
Study procedures were offered in English or Xhosa depending on the language preference of the mother.
This study was approved by the relevant university Health Research Ethics Committees
(University of Cape Town study number: 666/2021).
Informed consent was collected from mothers on behalf of themselves and their infants.
Demographic information, including maternal place of birth,
primary spoken language, maternal age at enrollment,
maternal educational attainment, and maternal income,
were collected at enrollment
(see @table1).

Families were invited to participate in three in-lab study visits
over their infant’s first two years of life.
At the first in-lab study visit (hereafter visit-1),
occurring when infants were between approximately 2 months and 6 months of age,
the following data were collected:
the infants' age (in months), sex,
infant electroencephalography (EEG),
and infant stool samples.

At the second study visit (hereafter visit-2),
occurring when infants were between approximately 6 months and 12 months of age
(age in months: M=8.60, SD=1.48, range=5.41-12.00)
and at the third study visit (hereafter visit-3),
occurring when infants were between approximately 12 months and 17 months of age
(age in months: M=14.10, SD=1.04, range=12.10-17.00),
infant EEG and stool samples were collected again.
At visits in which infants
were unable to complete both EEG and stool samples on the same day,
EEG and stool samples were collected on different days.
For concurrent time point analyses,
infants with EEG and stool collected more than two months apart were excluded.
Not all infants had EEG and microbiome data
collected at all three time points or contributed usable data at all three-time points.

All enrolled infants received a comprehensive medical exam at each visit,
which included assessments of eye-related conditions.
Several infants (n=3)
were identified as having eye-related anomalies during the medical exam,
and they were excluded from any further analyses.

=== EEG Processing

==== EEG Data Acquisition

Electroencephalography (EEG) data were acquired from infants while they were seated in their caregiver’s lap in a dimly-lit,
quiet room using a 128-channel high density HydroCel Geodesic Sensor Net
(EGI, Eugene, OR),
amplified with a NetAmps 400 high-input amplifier,
and recorded via an Electrical Geodesics, Inc.
(EGI, Eugene, OR)
system with a 1000 Hz sampling rate.
EEG data were online referenced to the vertex (channel Cz)
through the EGI Netstation software.
Impedances were kept below 100KΩ
in accordance with the impedance capabilities of the high-impedance amplifiers.
Geodesic Sensor Nets with modified tall pedestals
designed for improving inclusion of infants with thick/curly/tall hair were used as needed across participants
@mlanduEvaluatingNovelHighdensity2024.
Shea moisture leave-in castor oil conditioner
was applied to hair across the scalp prior to net placement
to improve both impedances and participant comfort
@mlanduEvaluatingNovelHighdensity2024.
This leave-in conditioner contains insulating ingredients,
so there is no risk of electrical bridging,
and has not been found to disrupt the EEG signal during testing (unpublished data).
Conditioning hair in this way allows for nets
to lay closer to the scalp for curly/coily hair types
and makes for more comfortable net removal at the end of testing.

The Visual-Evoked Potential (VEP) task
was presented using Eprime 3.0 software
(Psychology Software Tools, Pittsburgh, PA)
on a Lenovo desktop computer with an external monitor
19.5 inches on the diagonal facing the infant
(with monitor approximately 65 cm away from the infant).
A standard phase-reversal VEP was induced
with a black and white checkerboard
(1cm x 1 cm squares within the board)
stimulus that alternated presentation
(black squares became white, white squares became black)
every 500 milliseconds for a total of 100 trials.
Participant looking was monitored by video and by an assistant throughout data collection.
If the participant looked away during the VEP task,
the task was rerun.

==== EEG Data Pre-Processing

VEP data were exported from native Netstation `.mff` format
to `.raw` format and then pre-processed
using the HAPPE+ER pipeline within the HAPPE v3.3 software,
an automated open-source EEG processing software
validated for infant data
@monachinoHAPPEEventRelatedHAPPE2022.
A subset of the 128 channels
were selected for pre-processing that excluded the rim electrodes
as these are typically artifact-laden
(channels excluded from pre-processing included in Table S4).
The HAPPE pre-processing pipeline was run
with user-selected specifications outlined in Table S4.


Pre-processed VEP data were considered usable
and moved forward to VEP extraction if HAPPE pre-processing ran successfully,
at least 15 trials were retained following bad trial rejection,
and at least one good channel was kept within the visual ROI.
Note that channels marked bad during pre-processing
had their data interpolated as part of standard preprocessing pipelines for ERPs
@monachinoHAPPEEventRelatedHAPPE2022.
Interpolated channels were included in analyses here
as is typically done in developmental samples,
and given the low overall rates of interpolation present
(e.g., an average of between 4 to 5 of 5 possible good channels in the region of interest were retained at each visit time point).


==== Visual-Evoked Potentials (VEPs)

VEP waveforms were extracted and quantified
using the HAPPE+ER v3.3 GenerateERPs script
@monachinoHAPPEEventRelatedHAPPE2022.
Electrodes in the occipital region were selected as a region of interest
(i.e., E70, E71, E75, E76, E83).
The VEP waveform has three main components to be quantified:
a negative N1 peak,
a positive P1 peak,
and a negative N2 peak.
Due to normative maturation of the waveforms as infants age,
one set of user-specified windows for calculating component features
was used for visit-1 and 2 and another was used for visit-3.
For visits 1 and 2,
the window for calculating features for the N1 component was 40-100 ms,
75-175 ms for the P1 component,
and 100-325 ms for the N2 component.
For visit-3, the window for calculating features for the N1 component was 35-80 ms,
75-130 ms for the P1 component,
and 100-275 ms for the N2 component.
HAPPE+ER parameters used in extracting the ERPs are summarized in Table S5.

To correct for the potential influence of earlier components on later components,
corrected amplitudes and latencies were calculated and used in all analyses.
Specifically, the P1 amplitude was corrected for the N1 amplitude
(corrected P1 amplitude = P1 - N1 amplitude),
the P1 latency was corrected for the N1 latency
(corrected P1 latency = P1 - N1 latency),
the N2 amplitude was corrected for the P1 amplitude
(corrected N2 amplitude = N2 - P1 amplitude),
and the N2 latency was corrected for the P1 latency
(corrected N2 latency = N2 - P1 latency).

All VEPs were visually inspected to ensure that the automatically extracted values
were correct and were adjusted if observable peaks
occurred outside the automated window bounds.
Participants were considered to have failed this visual inspection
and were subsequently removed from the data set
if their VEP did not produce three discernible peaks.
VEP waveforms of included participants by time point are included in @figure2\A.
97 infants provided usable VEP data at visit-1,
130 infants provided usable VEP data at visit-2,
and 131 infants provided usable VEP data at visit-3.
For included participants,
EEG data quality metrics are summarized in Table S6.
T-tests for data quality metrics
(i.e., number of trials collected, number of trials retained,
number of channels retained in the ROI,
Pearson’s r for data pre- vs. post-wavelet thresholding
at 5, 8, 12, and 20 Hz)
were run between each visit combination
(i.e., visit-1 vs. visit-2, visit-1 vs. visit-3, visit-2 vs. visit-3).
For visits that differed in data quality,
follow-up post hoc correlations were run for the data quality measure
with each VEP feature at each visit in the T-test.
In no case did the data quality metric relate to VEP features at multiple visits,
making it highly unlikely the data quality difference contributed to results.

=== Biospecimens and sequencing

==== Sample Collection

Stool samples (n=315) were collected in the clinic
by the research assistant directly from the diaper,
transferred to Zymo DNA/RNA ShieldTM Fecal collection Tubes
(\#R1101, Zymo Research Corp., Irvine, USA)
and immediately frozen at -80 ˚C.
Stool samples were not collected
if the participant had taken antibiotics within the two weeks prior to sampling.

==== DNA Extraction

DNA extraction was performed at Medical Microbiology,
University of Cape Town, South Africa,
from stool samples collected in DNA/RNA Shield™ Fecal collection tube
using the Zymo Research Fecal DNA MiniPrep kit
(\# D4300, Zymo Research Corp., Irvine, USA)
following manufacturer’s protocol.
To assess the extraction process's quality,
ZymoBIOMICS® Microbial Community Standards
(\#D6300 and \#D6310, Zymo Research Corp., Irvine, USA)
were incorporated and subjected to the identical process as the stool samples.
The DNA yield and purity were determined using the NanoDrop® ND -1000
(Nanodrop Technologies Inc. Wilmington, USA).

==== Sequencing

Shotgun metagenomic sequencing was performed on all samples
at the Integrated Microbiome Research Resource
(IMR, Dalhousie University, NS, Canada).
A pooled library (max 96 samples per run)
was prepared using the Illumina Nextera Flex Kit
for MiSeq and NextSeq from 1 ng of each sample.
Samples were then pooled onto a plate
and sequenced on the Illumina NextSeq 2000 platform
using 150+150 bp paired-end P3 cells,
generating 24M million raw reads and 3.6 Gb of sequence per sample
@comeauPreparingMultiplexedWGS2023.

=== Statistics / computational analysis

==== Age-Related Changes in VEP Features

To determine age-related changes in VEP features,
six linear mixed models with each VEP feature as the outcome
(i.e., N1 amplitude/latency, P1 amplitude/latency, N2 amplitude/latency)
were run using the lme4 package
@batesFittingLinearMixedEffects2015
in R with age in months as the predictor of interest
and number of retained trials as a covariate.

==== Metagenome processing

Raw metagenomic sequence reads (2.5 x 10#super[7] ± 1.4 x 10#super[7] reads/sample)
were processed using tools from the bioBakery as previously described
@bonhamGutresidentMicroorganismsTheir2023
@beghiniIntegratingTaxonomicFunctional2021.
Briefly, KneadData v0.10.0 was used with default parameters
to trim low-quality reads and remove human sequences (using reference database hg37).
Next, MetaPhlAn v3.1.0 (using database mpa_v31_CHOCOPhlAn_201901)
was used with default parameters to map microbial marker genes to generate taxonomic profiles.
Taxonomic profiles and raw reads were passed
to HUMAnN v3.7 to generate stratified functional profiles.

==== Microbial community analysis

Principal coordinates analysis
was performed in the julia programming language
@bezansonJuliaFreshApproach2017
using the Microbiome.jl package
@bonhamMicrobiomeJlBiobakeryUtils2021.
Bray-Curtis dissimilarity (Distances.jl)
was calculated across all pairs of samples,
filtering for species-level classification.
Classical multidimensional scaling
was performed on the dissimilarity matrix (MultivariateStats.jl),
and axes with negative eigenvalues were discarded.

==== Feature Set Enrichment Analysis (FSEA)

Potentially neuroactive genesets
were extracted from Supplementary Dataset 1 from
@valles-colomerNeuroactivePotentialHuman2019.
Gut-brain modules provide Kegg Orthologue IDs (KOs)
@suzekUniRefComprehensiveNonredundant2007
@kanehisaKEGGResourceDeciphering2004,
which were mapped to UniRef90 IDs
using the utility mapping file provided with HUMAnN v3.1
@beghiniIntegratingTaxonomicFunctional2021.
For each stool/VEP pair,
logistic regression (LR) was performed,
linking the presence or absence of that UniRef in a sample
with each VEP feature (i.e., N1, P1, and N2 latencies and amplitudes),
controlling for the age at which the stool sample was collected,
the number of retained VEP trials,
and the difference in age between the stool collection and VEP measurement.
For concurrently collected stool and VEP comparisons (@figure2, @figureS2),
participants whose stool collection and VEP measurements
were more than two months apart were excluded.

$ "UniRef" ~ "vep" + "age_months" + "n_trials" + "age_diff" $

FSEA was performed on each geneset
that had at least five members in each comparison group
according to the procedure set out in Subramanian et. al. (2005)
@subramanianGeneSetEnrichment2005.
Briefly, enrichment scores (ES) are calculated
based on the rank order of z-statistics from the LR for each UniRef.
A permutation test was then performed
where the ES for 5000 random samples of ranks
of the same length as the gene set are calculated,
and the pseudo-p value is the fraction of permutations
where the permutation ES has a greater absolute value than the true ES.

Benjamini-Hochberg FDR correction
was performed separately on all concurrently tested geneset/VEP feature combinations
and all longitudinal geneset/VEP feature combinations.
Corrected p-values (q-values) less than 0.2
were considered statistically significant.

For longitudinal comparisons,
all participants that had a stool sample collected
at one visit and a VEP assessment at a subsequent visit were included
(visit-1 stool → visit-2 VEP, N = 84;
 v1 → v2, N = 76;
 v2 → v3, N = 69).
A total of 95 geneset/VEP features
were significant when using an FDR-corrected p-value cutoff of q < 0.2.
To ensure the robustness of these findings,
we randomly permuted participant IDs
between the stool and VEP assessments and repeated the analysis.
Over 10 random permutations,
a mean of 19.5 significant associations were identified,
suggesting that FDR correction is correctly calibrating the false-positive rate.

=== Data and code availability

Code for initial processing of data and for analyses performed in this manuscript
are available on github and archived on Zenodo
@bonhamKlepacCerajLabEEGMicrobiomeInitial2024.
Input data will be archived on Dryad
and downloadable via included scripts in analysis code upon publication.


== Results


=== The brain and microbiome develop rapidly in the first months of life

To investigate the co-development of the gut microbiome and visual neurodevelopment,
we collected stool and the VEP in a longitudinal cohort of 194 children in South Africa during the first 18 months of life
(@figure1\A, B, @table1;
visit-1, N = 119, age 3.6 ± 0.7 months,
visit-2, N = 144, age 8.7 ± 1.4 months,
visit-3, N = 130, age 14.2 ± 1.0 months).
As expected for children at this age,
both amplitude and latency VEP features were strongly correlated with age
@lippeElectrophysiologicalMarkersVisuocortical2007
@yassourNaturalHistoryInfant2016,
That is, as infants got older, N1 amplitude became more negative
(N1: b=-0.07, p<.05),
corrected P1 amplitude became smaller
(P1: b=-0.50, p<.05),
corrected N2 amplitude became smaller
(N2: b=0.52, p<.05),
and all latencies became shorter
(N1: b=-0.79, p<.05; P1: b=-1.21, p<.05; N2: b=-4.10, p<.05)
(@figure2\A).


Similarly,
microbial composition was developmentally dependent,
as expected
@koenigSuccessionMicrobialConsortia2011
@yassourNaturalHistoryInfant2016
@backhedDynamicsStabilizationHuman2015,
with early samples dominated by _Bifidobacterium_ and _Bacteroides_
species, while later samples have increasing _Prevotella_
and anaerobic genera such as _Faecalibacterium prausnitzii_
(@figure2\B).
Ordinations reveal a similar relationship with age,
with the first principal coordinate axis for both taxonomic profiles
(@figure2\C; variance explained = 15.1%; R = -0.50)
and functional profiles
(@figure2\D; variance explained = 12.9%; R = -0.57)
driven strongly by the age of the participant at the time of collection.

#figure(
  image("mainfigures/figure1.svg"),
  caption: [
    *Study design to capture the dynamic nature of early microbiome and brain development.*#linebreak()
    (A) Study design; participants (N=194) were seen up to 3 times over the first 18 months of life.
        Stool samples and EEG data were collected, generating microbial functional profiles (stool)
        and VEP waveforms (EEG) used in subsequent analyses.
    (B) Longitudinal sampling of study participants;
        Density plots (top) for stool and EEG collection show the ages represented in each visit.
        The scatter plot (bottom) shows individual participant visits.
        Dotted lines connect separate visits for the same participant.
        When stool and EEG data were collected for the same visit (purple) but not on the same day,
        dot represents the median age of collection,
        and vertical bars in blue and red represent stool and EEG collections respectively.
],
) <figure1>


#figure(
  image("mainfigures/figure2.svg"),
  caption: [
    *The gut microbiome and VEP both develop over the first 18 months of life.*#linebreak()
    (A) Individual VEP feature measurements for amplitude (left) and latency (right)
        for all participants and all visits in the study, separated by age. 
        mean and +/- standard error (S.E.) wave-forms for each visit are plotted below.
    (B) Relative abundance of the top 11 microbial species across all visits.
        All other species were summed so that the total abundance is 100%.
        each column represents a single sample, and samples are ordered by hierarchical clustering
        based on Bray-Curtis dissimilarity of the full microbial composition.
    (C) Principal coordinate analysis (PCoA) by multidimensional scaling (MDS)
        on Bray-Curtis dissimilarity of taxonomic profiles;
        percent variance explained (fraction of positive eigenvalues)
        by each of the first two axes are indicated on the x and y axes respectively.
    (D) PCoA of microbial functional profiles (UniRef90s).
],
) <figure2>


=== Microbial genes with neuroactive potential are associated with concurrently measured visual development

To test whether microbial metabolic potential was related to early life brain activity,
we performed feature set enrichment analysis (FSEA)
using previously-defined groups of potentially neuroactive microbial genes
and the concurrently measured VEP amplitude and latency features
@bonhamGutresidentMicroorganismsTheir2023
@valles-colomerNeuroactivePotentialHuman2019.
For each gene set that had at least 5 genes represented in a given comparison group,
logistic regression was performed using VEP features as predictors
and the presence or absence of each microbial gene in the metagenome
as the response to determine concurrent associations (see Methods).
Z statistics for in-set genes were compared to all genes
using a permutation test to determine significance of the associations
@subramanianGeneSetEnrichment2005.

#figure(
  image("mainfigures/figure3.svg"),
  caption: [
    *Feature Set Enrichment Analysis Reveals associations between microbial genes and VEP*#linebreak()
    (A) Volcano plots of of gene sets tested with feature set enrichment analysis (FSEA)
        for all 6 VEP features,
        with enrichment score (E.S.) compared to log scaled FDR-corrected p-value (Q).
        Colored dots were significantly enriched (positive E.S.) or depleted (negative E.S.)
        relative to the tested VEP feature.
    (B) Summary of results in (A), showing the fraction of each class of neuroactive genes
        (neurotransmitter metabolism, SCFA metabolism, amino acid metabolism, or other)
        that were statistically significantly enriched or depleted for each VEP feature
        for each visit in the analysis.
    (C) Enrichment plots for selected gene sets and there association with P1 latency.
        Each plot shows the distribution of associations of individual genes
        within the gene set and the VEP feature.
        Dots are colored if the geneset as a whole was significantly associated.
        Enrichment plots for all gene set / VEP feature associations may be found in Figure S2.
],
) <figure3>
#v(1em)

Of the 35 genesets assessed,
19 had sufficient representation to test,
and of those, 18 were significantly associated with at least one EEG feature
during at least one visit within the 18-month window,
after correcting for false discovery rate
(Benjamini-Hochberg, q < 0.2; @figure3, Table S2).
Microbial genes involved in synthesis or degradation of molecules with neuroactive potential across all categories considered
(i.e., neurotransmitters, amino acid metabolism, SCFAs, other)
were associated with both VEP amplitudes and latencies at each visit (Figure S2),
demonstrating widespread associations between early life gut microbiome and visual cortex neurodevelopment.
The number of these associations increased over time
(visit-1 had 6 associations, visit-2 had 24, and visit-3 had 37).

Specifically, across the gene sets involved in neurotransmitter synthesis and degradation,
glutamate synthesis/degradation and GABA synthesis showed associations with all VEP features,
primarily at the 2nd and 3rd visit
(mean ages 8.6 and 14.1 months, respectively; @figure3\C Table S2).
Gene sets involved in tryptophan metabolism and associated pathways
(i.e., quinolinic acid) were also strongly concurrently related to VEP development.
Specifically, tryptophan metabolism genes were associated with VEP latencies
just after each VEP component showed its greatest window of developmental change
(components emerge sequentially as follows: P1, N1, N2), that is,
of P1 at visit-1,
N1 and P1 at visit-2,
and N2 at visit-3.

Several short-chain fatty acid (SCFA)-metabolizing gene sets were also found to have multiple associations with VEP features.
Specifically, acetate synthesis was strongly associated with almost all VEP features
(Table S2).
Butyrate synthesis was associated with P1 and N2 amplitudes and latencies
around the end of the first year of life (visits 2 and 3),
when the visual cortex is most actively undergoing myelination.
Lastly, propionate synthesis/degradation was significantly associated with VEP latencies
at every visit over the 18-month window
(N1 at visits 1 and 2, and both P1 and N2 at visit-3).
These SCFA metabolizing genes showed almost double the concurrent associations
with VEP latencies than amplitudes
(11 associations with latencies, 6 associations with amplitudes).

Finally, within the remaining gene sets tested,
we observed robust associations in particular between menaquinone (Vitamin K2) gene sets
and the VEP features over this infancy window.
This is an expected relationship,
as vitamin K2 specifically is posited to promote healthy vision,
both outside of the brain through effects on the retina,
and within the brain where it can protect neural circuits from oxidative stress
@liNovelRoleVitamin2003.

Notably, across significant gene set associations with VEP features,
the P1 and N2 component amplitudes and latencies were consistently the most sensitive to these microbial gene sets.
Both P1 and N2 components are known to show the most dramatic changes with development during the first year of life
@lippeDifferentialMaturationBrain2009
and may best reflect underlying visual learning and plasticity at this stage
(Table S1).

=== Microbial metabolic potential predicts future brain development in infancy

We initially hypothesized that the earliest microbial influences
would have the largest effects on brain development,
but in cross-sectional analysis with concurrently measured VEP, we observed the fewest number of associations at visit-1. 
To differentiate whether this cross-sectional finding indicated the early microbiome
was sparsely related to visual cortical development at 4 months only or more broadly at all subsequent time points,
we sought to determine whether microbial genes at early time points were associated with subsequent VEP development.
To investigate this,
we performed FSEA on stool samples collected at visit-1 with visit-2 VEP
(see @table2, age at stool collection = 3.6 ± 0.8 months, age at VEP = 8.6 ± 1.5 months)
or visit-3 VEP
(see @table3, age at stool collection = 3.7 ± 0.7 months, age at VEP = 14.1 ± 1.1 months),
as well as visit-2 stool samples with visit-3 VEP
(see @table4, age at stool collection = 8.9 ± 1.5 months, age at VEP = 14.3 ± 1.0 months; Table S3)
(@figure4\A, Figure S3).


#figure(
  image("mainfigures/figure4.svg"),
  caption: [
    (A) Age distributions for cross-visit comparisons,
        with age at stool collection (left) and age of VEP measurement (right)
        for each participant included in the analysis.
        Collections for the same individual are connected by a dotted gray line.
    (B) Volcano plots of of gene sets tested with feature set enrichment analysis (FSEA)
        for all 6 VEP features,
        with enrichment score (E.S.) compared to log scaled FDR-corrected p-value (Q).
        Colored dots were significantly enriched (positive E.S.) or depleted (negative E.S.)
        relative to the tested VEP feature.
    (C) Summary of results in (B), showing the fraction of each class of neuroactive genes
        (neurotransmitter metabolism, SCFA metabolism, amino acid metabolism, or other)
        that were statistically significantly enriched or depleted for each VEP feature
        for each cross-visit comparison in the analysis.
    (D) Enrichment plots for selected gene sets and there association with P1 latency.
        Each plot shows the distribution of associations of individual genes
        within the gene set and the VEP feature.
        Dots are colored if the geneset as a whole was significantly associated.
        Enrichment plots for all gene set / VEP feature associations may be found in Figure S2.
],
) <figure4>
#v(1em)

All gene sets except those involved in the synthesis of 3,4-dihydroxyphenylacetic acid (DOPAC),
a metabolite of dopamine,
that had a significant hit with concurrently measured VEP
were also significantly associated with at least one future VEP feature
(@figureS3\B, Tables 2-5).
Notably, the quantity of those associations increased substantially
for all longitudinal comparisons compared to concurrent comparisons.
For example, only 6 visit-1 microbial gene sets were associated with visit-1 VEP,
and each of those was only associated with a single concurrently measured VEP feature.
By contrast, these longitudinal analyses revealed that visit-1 microbial gene sets
show a much richer pattern of associations with future VEP feature development.
Specifically, 13 visit-1 gene sets were associated with visit-2 VEP features,
and 11 were associated with visit-3 VEP features,
the majority
(9/13 for visit-2, 8/11 for visit-3)
were associated with at least 2 VEP features,
and nearly half
(6/13 for visit-2, 5/11 for visit-3)
were associated with more than 2 future VEP features.

Longitudinally, the early microbiome (visit-1)
was related to VEP features at visit-2 and visit-3 fairly evenly
(12/28 visit-1 microbiome associations to visit-3 VEP latencies,
 15/30 associations to visit-3 amplitudes),
suggesting early microbiome metabolism in the first 6 months of life is associated with visual neurodevelopment over the next year.
Microbiome metabolism from visit-2 was associated with
similar numbers of visit-3 VEP features as visit-1 microbiome
(17 visit-3 latency features, 18 visit-3 amplitude features),
suggesting continued co-development of these systems over the first postnatal year.
Neurotransmitters GABA and glutamate,
tryptophan metabolism (tryptophan and quinolinic acid),
SCFAs including acetate, butyrate, and propionate,
as well as menaquinone (Vitamin K2)
were again all significantly associated with multiple VEP features across multiple longitudinal comparisons.

Importantly, the nature and identity of the longitudinal associations
varied over development for many gene sets,
indicating temporal specificity to these associations.
With respect to the neurotransmitter-related pathways,
amongst the associations between GABA synthesis genes and future VEP features,
GABA genes specifically from visit-1 showed the greatest number of associations
with future VEP features (5/6 GABA associations) at visits 2 and 3 equally,
and the majority of these associations were VEP amplitudes
(reflecting development of neurotransmission including excitatory/inhibitory balance).
To a lesser extent, glutamate metabolism genes followed a similar temporal pattern
(8/13 glutamate associations involved visit-1 genes)
but did not relate to amplitudes or latencies differentially like GABA.
This pattern of results suggests early (within the first 6 postnatal months)
microbiome GABA/glutamate dynamics, especially GABA,
are most relevant for changes to visual cortex function over the following year.

Tryptophan-related pathway genes (responsible for generating serotonin,
amongst other products) from visit-1 were also responsible
for the majority of associations with future VEP features
(tryptophan: 6/8 associations; quinolinic acid: 6/9 associations).
In contrast to GABA but similar to glutamate,
tryptophan-related gene associations were largely shorter-term associations
with VEP features at the visit immediately following gene set measurement
(approximately 5 months later; 12/17 associations),
indicating dynamic co-development over the first 18 months of life.
Across neurotransmitter-related gene set associations
(GABA, glutamate, tryptophan/serotonin),
there was thus a clear pattern whereby early (\~4 months old)
microbiome gene sets showed the largest number of associations with subsequent VEP feature development.

SCFAs showed a different developmental pattern of associations with future VEP features.
Specifically, propionate and butyrate metabolism genes
from both visit-1 and visit-2 showed associations with future VEP features,
but here the effects were almost entirely observed
for VEP features at visit-3 (10/10 propionate and 7/8 butyrate associations).
Moreover, acetate and butyrate metabolism genes
were doubly associated with future VEP latencies compared to amplitude features, possibly reflecting their known association with myelination pathways.

Finally, menaquinone (Vitamin K2) metabolism genes
followed a similar pattern to the SCFAs
in that genes from visit-1 and visit-2 were largely associated
with future VEP features at visit-3 (8/9 menaquinone associations).
This indicates persistent associations of this early microbiome gene set
with individual differences in VEP features early in the second year of life.
It is possible that extrinsic factors related to development
mutually influence both the gut microbiome and neural development,
though we additionally tested whether VEP features
were associated with microbial metabolism at a future visit,
and found substantially fewer associations (29 total associations,
compared with 95 when analyzing early stool samples with future VEP).
While this does not prove a causal relationship,
it is consistent with the hypothesis that microbial metabolism influences brain development.

== Discussion

The past decade has seen a remarkable growth in our understanding of
the relationships between the developmental changes of the gut microbiome and the brain.
However, a great deal of that investigation has focused on adult populations
or neuropsychiatric disorders, limiting the potential to explain
how these associations emerge during early development.
Here, we address this key open question
by leveraging a rich longitudinal dataset over the first year of life,
which is the time of greatest developmental change
for both the microbiome and brain given the unfolding of foundational sensory neurodevelopment.
Our data revealed that microbial genes involved in the metabolism of neuroactive molecules
are associated with concurrent and subsequent visual cortical neurodevelopment.
These pathways included those for the neurotransmitters GABA and glutamate,
the amino acid tryptophan,
and short-chain fatty acids involved in myelination,
including acetate and butyrate.


Specifically, we have shown a robust, prospective relationship
between microbial genes involved in the metabolism of neuroactive compounds
and the development of visual cortical function as measured by the VEP electrophysiological response.
We found that microbial metabolism
is more strongly associated with future measures of the VEP
than those collected concurrently.
While not dispositive, this would be the predicted outcome
if microbial genes are causally influencing brain development.
As additional evidence for this interpretation,
we did not observe the same rich set of associations
in the converse analyses examining whether VEP related to future microbiome properties.
Microbial metabolism within the first 6 months
shows the most associations with subsequent visual neurodevelopment,
suggesting the early postnatal microbiome may play a particularly important role
in the co-development of these systems.
This interpretation is also supported by prior research
showing that associations of the microbiome with behavioral readouts of neurocognition
are stronger prospectively than concurrently
@carlsonInfantGutMicrobiome2018b.
Moreover, specific associations between gene sets and VEP features
showed temporal specificity within the 18-month developmental window assessed,
suggesting that the impact of early microbial metabolism on the brain is developmentally dependent.


Notably, the gene sets most highly associated
with visual functional neurodevelopment over infancy
are for the metabolism of molecules with known links to developmental neuroplasticity
@fagioliniInhibitoryThresholdCriticalperiod2000
@henschLocalGABACircuit1998
@takesianInhibitoryCircuitGating2018.
Specifically, we observed associations for gene sets related to glutamate and GABA,
neurotransmitters that are central to regulating excitatory/inhibitory (E/I) cortical balance.
Developmental changes in E/I balance modulate the degree of neuroplasticity in the mammalian cortex,
including regulating the start and progression
of critical period neuroplasticity mechanisms in the visual cortex
@henschLocalGABACircuit1998
@fagioliniInhibitoryThresholdCriticalperiod2000
@gabard-durnamSensitivePeriodsHuman2020
@margolisLongitudinalEffectsPrenatal2024.
Our observed pattern of results suggests early
(within the first 6 postnatal months)
microbiome GABA/glutamate dynamics, especially GABA,
are most relevant for changes to visual cortex function over the following year.
Gut production of GABA may influence cortical GABA levels
via active transport from bloodstream to brain
@takanagaGAT2BGT1System2001
@al-sarrafTransport14CgaminobutyricAcid2002
@shyamaladeviEvidenceThatNitric2002.
Recent evidence suggests gut-derived glutamate
may also influence brain levels and function
@matsumotoCerebralLowmolecularMetabolites2013
@kawaseGutMicrobiotaMice2017
@filpaRoleGlutamatergicNeurotransmission2016
and can operate via indirect mechanisms
(either transformation into GABA
or via regulating glutamate levels in the bloodstream
that impact glutamate transfer from brain to bloodstream).


Tryptophan related pathway genes were also identified here
that are responsible for generating serotonin
as well as other neuroactive molecules such as kynurenic acid (an SMDAR antagonist)
@agusGutMicrobiotaRegulation2018.
Both serotonin and kynurenic acid
are implicated in early neuroplasticity and neurotransmitter regulation,
and serotonin has potent effects for visual cortex plasticity in particular
@vetencourtSerotoninTriggersTransient2011
@guInvolvementSerotoninDevelopmental1995.
While quinolinic acid is part of the kynurenine pathway
and is a neurotoxin that can cause neuronal dysfunction,
it may also play a role in glutamate uptake in the brain
@lugo-huitronQuinolinicAcidEndogenous2013
@tavaresQuinolinicAcidInhibits2000.
Specifically, tryptophan metabolism genes
were associated with VEP latencies just after each VEP component showed its greatest window of developmental change
(components emerge sequentially as follows: P1, N1, N2),
that is, of P1 at visit-1, N1 and P1 at visit-2, and N2 at visit-3.
This pathway may thus relate to processes stabilizing the neural circuitry
(i.e. downregulating neuroplasticity)
underlying each VEP component,
an account consistent with recently observed effects
of serotonin within the visual cortex in rodents
@hongNorepinephrinePotentiatesSerotonin2022.
Importantly, nearly all of the body’s serotonin
is produced in the gut by enterochromaffin cells,
and this biosynthesis is regulated by microbes
@yanoIndigenousBacteriaGut2015
@jamesonLinkingGutMicrobiota2018,
making this pathway an especially promising candidate intervention target
for further research in development.


We further found that gene sets for short-chain fatty acids
important in downregulating neuroinflammation and promoting myelination within the brain
were robustly related to visual neurodevelopment.
Myelination is important for down-regulating plasticity in neural circuitry over development
by stabilizing and protecting circuits that have been shaped by early experience
@takesianBalancingPlasticityStability2013.
Specifically, we observed associations between
acetate, butyrate, and proprionate genes with VEP development.
Acetate is a critical component required for the increased lipid synthesis
that happens during postnatal myelination in the brain
@madhavaraoDefectiveNacetylaspartateCatabolism2005.
Circulating butyrate also increases myelination
@chenButyrateSuppressesDemyelination2019,
and though propionate’s relation to myelinating oligodendrocytes remains unclear,
it is known to protect myelinating Schwann cells outside of the brain
from oxidative stress
@gruterPropionateExertsNeuroprotective2023.
Acetate, butyrate, and propionate are all also widely regarded as neuroprotective
by promoting healthy microglial development and downregulating neuroinflammation
that interferes with myelination
@caetano-silvaInhibitionInflammatoryMicroglia2023.
VEP latency features reflect myelination
@margolisLongitudinalEffectsPrenatal2024
@barnikolPatternReversalVisual2006
@youLatencyDelayVisual2011,
and accordingly, these SCFAs showed more associations
with VEP latency features prospectively,
especially VEP latency features in vsit 3.
This pattern of results is consistent with these SCFA roles’ in myelination
occurring over the second half of the developmental window studied.
SCFAs including acetate, butyrate, and propionate
can pass the blood-brain barrier to directly influence myelination-related processes
within the brain.
Taken together, the pattern of results
across GABA, glutamate, tryptophan, and SCFA gene sets
suggests early postnatal microbiome-derived metabolites
relate to key neuroplasticity regulation processes within the cortex.

Our study is a substantial advancement over prior work
on the microbial-gut-brain axis in early life
due to the sequencing method,
the large number of participants,
the longitudinal study design,
and the inclusion of participants from a previously under-represented region of the world.
The use of shotgun metagenomic sequencing
enables direct interrogation of microbial metabolic potential.
Prior research primarily used amplicon (16S rRNA gene) sequencing,
which enables lower-resolution taxonomic identification
and is restricted to inferring metabolic potential based on taxonomy.
Moreover, while several studies in infancy have inferred gut-brain associations
by linking microbiome measures to subsequent neurodevelopmental measures using behavioral assessments
(e.g., Bayley Scales of Infant Development),
noting associations with visually-mediated cognition
@carlsonInfantGutMicrobiome2018b,
this study assessed the development of the visual cortex directly using the VEP.
The VEP is advantageous
because features reflect largely neurotransmission-related (via amplitudes)
or structural (i.e., myelination, via latencies)
changes over this developmental window, facilitating some specificity
in the observed associations.
Moreover, the VEP can be indexed with fidelity from birth,
providing a continuous measure of visual cortical function across the study age-range.
Additionally, this study involved a large number of participants (194)
contributing dense longitudinal data,
with up to three time points,
all taken in the first 18 months of an infant’s life.
While prior work focused on single time point measures of microbiome and neurodevelopment
@gaoGutMicrobiomeBrain2019,
longitudinal associations allowed us to investigate the changing relation between gut microbial metabolism and the development of visual neurocircuitry over time.

One limitation of this study is the fact that we are only able to observe the genomic composition of the microbiome,
rather than the concentration of metabolites themselves.
This prevents us from determining the concentration of these molecules
in the gastrointestinal tract, blood, and brain,
as the abundance of these genes does not provide information about their activity,
their interactions with other metabolic pathways (including those of the host),
or absorption by colonic epithelial cells.
Moreover, the relationship between gene abundance and molecule concentration
may be counterintuitive, since the relationship between degradation and synthesis
of metabolites occurs both at the individual organism level and at the community level.
For example, genes for breaking down a molecule may be prevalent if that molecule is at high concentrations,
or the molecules may be rapidly degraded
by other members of the community the moment they are produced.
Furthermore, it may be that the relation between metabolite and brain development remains stable over time,
but the relation between molecules and microbial
selection changes at different stages of life.
Addressing these limitations in humans is challenging,
even if looking at stool metabolites,
because overall exposure throughout the gastrointestinal tract
is not necessarily reflected in the final concentration of those molecules in stool.
Therefore, metabolites from blood plasma
could provide more accurate systemic concentrations of molecules,
but challenges remain on how to interpret them in humans
@dengComparisonFecalBlood2023
@dekkersOnlineAtlasHuman2022.

Given that the VEP is evolutionarily conserved in mammals
and can be accurately measured during development,
the hypotheses generated in humans in this study
are readily testable mechanistically using in vivo models in future research.
For example, VEP could be assessed in germ-free or defined-microbiome animals
(Wymore Brand et al. 2015; Kennedy, King, and Baldridge 2018)
and may be supplemented with specific molecules such as SCFAs,
or colonized with microbial species lacking or providing specific metabolic pathways.
Furthermore, molecule concentrations in tissues
from the gut to the brain can be directly assessed in these models.
Uncovering relations between microbial metabolism and specific molecules
may also generate hypotheses that can be confirmed in human data.
This study, therefore,
provides a foundation for deep investigation
of the link between the human gut microbiome and brain development.

== Acknowledgments

We would like to extend our thanks
to all the caregivers and their infants
who generously provided their time and samples for this study.
We are also grateful to the dedicated nurses and researchers
who recruited participants and collected data,
ensuring the success of this project.
The Khula South African Data Collection Team is
Layla Bradford, Simone Williams, Lauren Davel,
Tembeka Mhlakwaphalwa, Bokang Methola, Khanyisa Nkubungu,
Candice Knipe, Zamazimba Madi, and Nwabisa Mlandu.
This research was supported by the Wellcome LEAP 1kD program.

=== CRediT

Conceptualization: KSB, ETM, KAD, LJGD, VKC#linebreak()
Data curation: KSB, ETM, GFB, FP, SM, MRZ, MM, DH#linebreak()
Formal analysis: KSB, ETM#linebreak()
Funding acquisition: KSB, NP, DCA, DKJ, SCRW, DA, MG, WPF, KAD, LJGD, VKC#linebreak()
Investigation: KSB, ETM, GFB, SM, MRZ, MM, DH#linebreak()
Methodology: KSB, NP, DCA, DKJ, SCRW, DA, MG, WPF, KAD, LJGD, VKC#linebreak()
Project administration: MRZ, MM, DH, LJGD, VKC#linebreak()
Software: KSB, ETM, GFB, CH, LJGD#linebreak()
Resources: SM, AS, FP, MRZ, MM, LD, CB, Khula, KAD, VKC#linebreak()
Supervision: KSB, LJGD, VKC#linebreak()
Validation: KSB, ETM#linebreak()
Visualization: KSB, ETM#linebreak()
Writing – original draft: KSB, ETM, LJGD, VKC#linebreak()
Writing – review & editing: KSB, ETM, CH, NP, KAD#linebreak()

== Disclosures

The authors declare no financial or other conflicts of interest.

== Author correspondence

#grid( columns: (1fr, 1fr), align: center, gutter: 4pt,
      [ Vanja Klepac-Ceraj ],      [ Laurel Gabard-Durnam ],
      [ vklepacc\@wellesley.edu ], [ l.gabard-durnam\@northeastern.edu ],
      [ 781-283-3541 ],            [ 628 ISEC ],
      [ 106 Central St ],          [ Northeastern University ],
      [ Wellesley, MA 02481] ,     [Boston, MA 02482 ]
)



#bibliography("refs.bib", style: "biological-psychiatry.csl")

#show figure.where(
  kind: table
): set figure.caption(position: top)

#show figure.caption: it => context [
  *#it.supplement~#it.counter.display()#it.separator*#it.body
]

#pagebreak()

== Tables


#set text(8pt)
#show table.cell: it => {
  if it.x == 1 {
    set align(center)
    it
  } else {
    set align(left)
    it
  }
}
#figure(
    table(columns: (4fr, 1fr),
          inset: 3pt,
          stroke: none,
            table.hline(stroke: 2pt),
            [], [*Overall\
                (N=194)*],

            table.hline(),
            [*Mean (SD) Age at EEG Data Collection (months)*], [],
            [Visit 1 (N=97)],     [3.7 (0.85)],
            [Visit 2 (N=129)],           [8.6 (1.46)],
            [Visit 3 (N=130)],           [14.1 (1.03)],

            [*Mean (SD) Age at Stool Data Collection (months)*], [],
            [Visit 1 (N=119)],           [3.6 (0.76)],
            [Visit 2 (N=105)],           [8.8 (1.43)],
            [Visit 3 (N=91)],            [14.0 (1.24)],

            [*Maternal Place of Birth*], [],
            [South Africa],              [191 (98.5%)],
            [In the African Continent (not South Africa)],  [3 (1.5%)],

            [*Primary Spoken Language*], [],
            [Xhosa Language],            [187 (96.4%)],
            [Sotho Language],            [2 (1.0%)],
            [Zulu Language],             [1 (0.5%)],
            [English Language],          [2 (1.0%)],
            [Ndebele Language],          [1 (0.5%)],
            [Afrikaans Language],        [1 (0.5%)],

            [*Maternal Age at Infant Birth (years)*], [],
            [Mean (SD)],                 [29.2 (5.63)],
            [Median [Min, Max]],         [29.0 [18.0, 41.0]],
            [Missing],                   [1 (0.5%)],

            [*Maternal Educational Attainmentᵃ*], [],
            [Completed Grade 6 (Standard 4)
             to Grade 7 (Standard 5)],           [4 (2.1%)],
            [Completed Grade 8 (Standard 6)
             to Grade 11 (Standard 9) i.e.,
             high school without matriculating], [78 (40.2%)],
            [Completed Grade 12 (Standard 10)
             i.e., high school]       ,          [88 (45.4%)],
            [Part of university/ college/
             post-matric education],             [13 (6.7%)],
            [Completed university/ college/
             post-matric education],             [11 (5.7%)],

            [*Maternal Monthly Incomeᵇ (South African Rand/ZAR)*], [],
            [Less than R1000 per month],          [97 (50.0%)],
            [R1000 - R5000 per month],            [76 (39.2%)],
            [R5000 - R10,000 per month],          [16 (8.2%)],
            [More than R10,000 per month],        [0 (0%)],
            [Unknown],                            [5 (2.6%)],

            [*Infant Biological Sex*], [],
            [Female],            [91 (46.9%)],
            [Male],              [103 (53.1%)],
            table.hline(stroke: 2pt),

            table.cell(colspan:2)[ᵃThe South African Educational System was formerly divided into years called standards,
                    similarly to the way the United States Educational System is divided into grades.
                    The equivalent in terms of standards is provided in parentheses next to each mentioned grade.
                    “University/College/Post-Matric Education” refers to tertiary
                    or post-secondary education as defined by the World Bank.],
            table.cell(colspan:2)[ᵇAt the time of writing (1/16/24), 1 US Dollar = 18.87 South African Rand (ZAR).]
            ),
    caption: [*Overall Demographic Information*]
) <table1>


#set table.cell(align: center+horizon)


#figure(
    caption:[*Longitudinal FSEA, visit-1 stool #sym.arrow.r.stroked visit-2 VEP*],
    table(
        columns: 5,
        stroke: 0.5pt,
        inset: 3pt,
        table.hline(stroke: 2pt),
    [Gene set], [Feature], [Component], [Enrichment], [Q value],
            table.hline(),
    table.cell(rowspan:2, /*fill: rgb("#89CDD8")*/)[GABA synthesis],
        table.cell(rowspan:2)[amplitude],
            [P1], [0.3987], [0.0807],
            [N2], [0.3661],  [0.1129],
    table.cell(rowspan:3, /*fill: rgb("#89CDD8")*/)[Glutamate synthesis],
        table.cell(rowspan:2)[latency],
            [N1], [0.2829], [0.0445],
            [P1], [0.2625], [0.0807],
        [amplitude],
            [P1], [0.2316], [0.1023],
    table.cell(rowspan:3, /*fill: rgb("#89CDD8")*/)[Glutamate degradation],
        table.cell(rowspan:2)[latency],
            [P1], [0.3369], [0.1598],
            [N2], [0.3252], [0.1800],
        [amplitude],
            [P1], [0.3256], [0.1760],
    table.cell(rowspan:4, /*fill: rgb("#82B574")*/)[Tryptophan synthesis],
        table.cell(rowspan:3)[latency],
            [N1], [0.2211], [0.0272],
            [P1], [0.1671], [0.0807],
            [N2], [0.1900], [0.0484],
        [amplitude],
            [N2], [0.1630], [0.1023],
    table.cell(/*fill: rgb("#82B574")*/)[Quinolinic acid synthesis],
        [amplitude],
            [N1], [0.2735], [0.1609],
    table.cell(rowspan:2, /*fill: rgb("#82B574")*/)[Quinolinic acid degradation],
        table.cell(rowspan:2)[latency],
            [P1], [0.1979], [0.1023],
            [N2], [0.1950], [0.1118],
    table.cell(rowspan:5, /*fill: rgb("#F29972")*/)[Acetate synthesis],
        table.cell(rowspan:3)[latency],
            [N1], [0.2190], [0.0293],
            [P1], [0.1926], [0.0611],
            [N2], [0.2157], [0.0318],
        table.cell(rowspan:2)[amplitude],
            [P1], [0.1860], [0.0742],
            [N2], [0.1960], [0.0549],
    table.cell(/*fill: rgb("#F29972")*/)[Butyrate synthesis],
        [amplitude],
            [N1], [0.3193], [0.1698],
    table.cell(rowspan:3, /*fill: rgb("#F29972")*/)[Isovaleric acid synthesis],
        [latency],
            [N1], [0.2711], [0.1223],
        table.cell(rowspan:2)[amplitude],
            [N1], [0.2412], [0.1800],
            [P1], [0.2409], [0.1609],
    table.cell(/*fill: rgb("#D0D17D")*/)[Menaquinone synthesis],
        [latency],
            [N1], [0.1393], [0.1598],
    table.cell(rowspan:2, /*fill: rgb("#D0D17D")*/)[Inositol synthesis],
        [amplitude],
            [N1], [0.4295], [0.1804],
        [latency],
            [N1], [0.4219], [0.1834],
    table.cell(rowspan:3, /*fill: rgb("#D0D17D")*/)[p-Cresol synthesis],
        table.cell(rowspan:3)[amplitude],
            [N1], [0.3665], [0.1825],
            [P1], [0.4520], [0.1014],
            [N2], [0.5164], [0.0445],
    table.cell(/*fill: rgb("#D0D17D")*/)[17-beta-Estradiol degradation],
        [latency],
            [N2], [0.3550], [0.1800],
    table.hline(stroke: 2pt),
    )
) <table2>

#figure(
    caption:[*Longitudinal FSEA, visit-1 stool #sym.arrow.r.stroked visit-3 VEP*],
    table(
        columns: 5,
        stroke: 0.5pt,
        inset: 3pt,
        table.hline(stroke: 2pt),
    [Gene set], [Feature], [Component], [Enrichment], [Q value],
        table.hline(),

    table.cell(rowspan:3, /*fill: rgb("#89CDD8")*/)[GABA synthesis],
        [latency],
            [P1], [0.4343], [0.0445],
        table.cell(rowspan:2)[amplitude],
            [P1], [0.4099], [0.0742],
            [N2], [0.5065], [0.0254],
    table.cell(rowspan:2, /*fill: rgb("#89CDD8")*/)[Glutamate synthesis],
        [latency],
            [N1], [0.3035], [0.0293],
        [amplitude],
            [N1], [0.2344], [0.1118],
    table.cell(rowspan:2, /*fill: rgb("#82B574")*/)[Tryptophan synthesis],
        [latency],
            [N1], [0.1875], [0.0445],
        [amplitude],
            [P1], [0.2176], [0.0159],
    table.cell(rowspan:3, /*fill: rgb("#82B574")*/)[Quinolinic acid degradation],
        table.cell(rowspan: 2)[latency],
            [P1], [0.1862], [0.1093],
            [N2], [0.2624], [0.1834],
        [amplitude],
            [P1], [0.1689], [0.1609],
    table.cell(/*fill: rgb("#F29972")*/)[Acetate synthesis],
        [latency],
            [N1], [0.1728], [0.0807],
    table.cell(/*fill: rgb("#F29972")*/)[Propionate synthesis],
        [amplitude],
            [P1], [0.3586], [0.0807],
    table.cell(rowspan: 4, /*fill: rgb("#F29972")*/)[Propionate degradation],
        table.cell(rowspan: 2)[latency],
            [N1], [0.5097], [0.1093],
            [P1], [0.5345], [0.0880],
        table.cell(rowspan: 2)[amplitude],
            [P1], [0.7388], [0],
            [N2], [0.6669], [0.0272],
    table.cell(rowspan: 3, /*fill: rgb("#F29972")*/)[Butyrate synthesis],
        table.cell(rowspan: 2)[latency],
            [N1], [0.3373], [0.1223],
            [N2], [0.3232], [0.1498],
        [amplitude],
            [P1], [0.3688], [0.1014],
    table.cell(rowspan:4, /*fill: rgb("#D0D17D")*/)[Menaquinone synthesis],
        [latency],
            [P1], [0.1582], [0.0870],
        table.cell(rowspan: 3)[amplitude],
            [N1], [0.1417], [0.1397],
            [P1], [0.1488], [0.1131],
            [N2], [0.1567], [0.1023],
    table.cell(/*fill: rgb("#D0D17D")*/)[Inositol synthesis],
        [amplitude],
            [N2], [0.4538], [0.1210],
    table.cell(rowspan: 3, /*fill: rgb("#D0D17D")*/)[ClpB],
        [latency],
            [N1], [0.2835], [0.1804],
        table.cell(rowspan: 2)[amplitude],
            [P1], [0.3284], [0.1129],
            [N2], [0.3585], [0.0800],
    table.hline(stroke: 2pt),
    )
) <table3>

#figure(
    caption:[*Longitudinal FSEA, visit-2 stool #sym.arrow.r.stroked visit-3 VEP*],
    table(
        columns: 5,
        stroke: 0.5pt,
        inset: 3pt,
        table.hline(stroke: 2pt),
    [Gene set], [Feature], [Component], [Enrichment], [Q value],
    table.cell(/*fill: rgb("#89CDD8")*/)[GABA synthesis],
        [latency],
            [N1], [0.3762], [0.1073],
    table.cell(rowspan:3, /*fill: rgb("#89CDD8")*/)[Glutamate synthesis],
        table.cell(rowspan:2)[latency],
            [N1], [0.1676], [0.1609],
            [P1], [0.2979], [0.0159],
        [amplitude],
            [P1], [0.2342], [0.0381],
    table.cell(rowspan:2, /*fill: rgb("#89CDD8")*/)[Glutamate degradation],
        [latency],
            [P1], [0.3174], [0.1404],
        [amplitude],
            [N2], [0.2829], [0.1875],
    table.cell(rowspan:2, /*fill: rgb("#82B574")*/)[Tryptophan synthesis],
        table.cell(rowspan:2)[amplitude],
            [N1], [0.1554], [0.0608],
            [P1], [0.1865], [0.0293],
    table.cell(rowspan:3, /*fill: rgb("#82B574")*/)[Quinolinic acid synthesis],
        table.cell(rowspan:2)[latency],
            [N1], [0.3186], [0.0610],
            [P1], [0.3331], [0.0437],
        [amplitude],
            [P1], [0.2577], [0.1391],
    table.cell(rowspan: 4, /*fill: rgb("#F29972")*/)[Acetate synthesis],
         table.cell(rowspan:3)[latency],
            [N1], [0.1673], [0.0610],
            [P1], [0.2574], [0.000],
            [N2], [0.1373], [0.1529],
        [amplitude],
            [P1], [0.1469], [0.1210],
    table.cell(/*fill: rgb("#F29972")*/)[Propionate synthesis],
        [amplitude],
            [P1], [0.2994], [0.1262],
    table.cell(rowspan: 4, /*fill: rgb("#F29972")*/)[Propionate degradation],
        table.cell(rowspan:2)[latency],
            [N1], [0.5503], [0.0742],
            [N2], [0.4492], [0.1609],
        table.cell(rowspan:2)[amplitude],
            [P1], [0.5088], [0.1118],
            [N2], [0.4229], [0.1875],
    table.cell(rowspan: 4, /*fill: rgb("#F29972")*/)[Butyrate synthesis],
        table.cell(rowspan:3)[latency],
            [N1], [0.3316], [0.0807],
            [P1], [0.3335], [0.0742],
            [N2], [0.2908], [0.1223],
        [amplitude],
            [N2], [0.3066], [0.1118],
    table.cell(rowspan: 2, /*fill: rgb("#F29972")*/)[Isovaleric acid synthesis],
        [latency],
            [P1], [0.2316], [0.1556],
        [amplitude],
            [P1], [0.2469], [0.1129],
    table.cell(rowspan: 4, /*fill: rgb("#D0D17D")*/)[Menaquinone synthesis],
        [latency],
            [N1], [0.1528], [0.1210],
        table.cell(rowspan:3)[amplitude],
            [N1], [0.1644], [0.1073],
            [P1], [0.2014], [0.0293],
            [N2], [0.1634], [0.1121],
    table.cell(rowspan: 3, /*fill: rgb("#D0D17D")*/)[Inositol degradation],
        table.cell(rowspan:2)[latency],
            [N1], [0.6381], [0.0293],
            [P1], [0.5323], [0.1014],
        [amplitude],
            [P1], [0.6413], [0.0293],
    table.cell(rowspan: 3, /*fill: rgb("#D0D17D")*/)[p-Cresol synthesis],
        table.cell(rowspan:2)[latency],
            [P1], [0.3118], [0.1397],
            [N2], [0.3152], [0.1391],
        [amplitude],
            [N1], [0.3039], [0.1658],
    table.cell(rowspan: 2, /*fill: rgb("#D0D17D")*/)[S-Adenosylmethionine synthesis],
        table.cell(rowspan:2)[amplitude],
            [N1], [0.2228], [0.1556],
            [P1], [0.2295], [0.1347],
    table.cell(/*fill: rgb("#D0D17D")*/)[17-beta-Estradiol degradation],
        [amplitude],
            [P1], [0.2818], [0.1875],
    table.cell(/*fill: rgb("#D0D17D")*/)[ClpB],
        [latency],
            [P1], [0.2687], [0.1804],
    )
) <table4>

#pagebreak()

== Supplemental Figures

#figure(
  image("mainfigures/figureS2.svg"),
  caption: [
    *Concurrent feature set enrichment analysis of microbial neuroactive genes and VEP for three visits.*
    FSEA results for all genesets where at least one visit had a significant hit (q < 0.2)
    with at least one VEP latency (A) or amplitude (B).
    Dots indicate the Z-statistic from logistic regression for each gene in a gene set.
    Vertical bars indicate the median Z-statistic for the gene set as a whole.
    Y-axis position for each gene set indicates visit number.
    Visit 1 for inositol degradation and DOPAC synthesis were not tested,
    since there were fewer than 5 genes from those genesets present in the sample (See Methods).
],
) <figureS2>

#figure(
        image("mainfigures/figureS3.svg"),
        caption: [
        *Gut microbial genes predict future VEP latencies and amplitudes.*
        (A) Age distributions for stool samples (left) and VEP (right)
        for each longitudinal comparison (same individual) tested,
        V1 stool →V2 VEP, V1 stool →V3 VEP, and V2 stool →V3 VEP.
        As in Figure S2, (B) and (C) show FSEA results for all genesets
        where at least one visit had a significant hit (q < 0.2)
        with at least one VEP latency or amplitude respectively.
        Dots indicate the Z-statistic from logistic regression for each gene in a gene set.
        Vertical bars indicate the median Z-statistic for the gene set as a whole.
        The Y-axis position for each gene set indicates longitudinal comparison.
        V1 → V2 and V1 → 3 for inositol degradation and DOPAC synthesis were not tested,
        since there were fewer than 5 genes from those genesets present in the sample (See Methods).
]) <figureS3>


