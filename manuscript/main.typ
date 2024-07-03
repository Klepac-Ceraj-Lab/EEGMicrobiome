#import "template.typ": *

// Take a look at the file `template.typ` in the file panel
// to customize this template and discover how it works.
#show: project.with(
  title: "Early infant microbiobial metabolism alters brain development",
  authors: (
    "Kevin S. Bonham",
    "Emma Margolis",
    "Laurel Gabard-Durnam",
    "Vanja Klepac-Ceraj",
  ),
  abstract: [
    Infancy is a time of rapid brain development when foundational neurologically-embedded sensory learning gut microbiome,
    also undergoing extensive developmental changes in early life,
    may alter brain development through the metabolism of neuroactive compounds.
    Here, we show across the first 18 months of life,
    microbial genes encoding enzymes that produce and degrade neuroactive compounds,
    including short-chain fatty acids and neurotransmitters GABA and glutamate,
    are associated with visual neurodevelopmental learning,
    measured by the visual-evoked potential (VEP).
    Concurrently measured microbial gene sets were most associated with VEP features when children were approximately 14 months old,
    but genesets from stool collected around 4 months of age were strongly associated with VEP features measured at 9 or 14 months of age,
    suggesting microbial metabolism in early life may have long term effects on neural development in a major sensory domain.
   ],
)

// Leave comments like this

= Introduction

The gut microbiome in early life has potential long-term implications for health and neurodevelopment.
One important way this influence can occur is through interactions with the central nervous system as a “microbial-gut-brain-axis”
@rheePrinciplesClinicalImplications2009
(Rhee, Pothoulakis, and Mayer 2009).
The gene-functional potential of the microorganisms that inhabit the gut exceeds that of human genomes by a hundredfold
and includes the ability to metabolize and synthesize many potentially neuroactive compounds
(Valles-Colomer et al. 2019).
Extensive work in preclinical models suggests that these potentially neuroactive compounds
can influence the brain through both direct and indirect routes.
For example, major neurotransmitters (e.g. glutamate, γ-aminobutyric acid (GABA), serotonin, dopamine)
are readily synthesized and degraded by intestinal microbes
and can directly stimulate enteric nerve cells or enter circulation and pass the blood-brain barrier
(Janik et al. 2016).
Glutamatergic/GABA-ergic signaling is critical for balancing the brain’s excitatory and inhibitory neurotransmission levels,
and bi-directional glutamatergic/GABA-ergic signaling between gut microbiome and brain
is implicated in several physical and mental health conditions
(Filpa et al. 2016; Baj et al. 2019).
Similarly the gut microbiome is critical to neurotransmitter metabolism regulation for serotonin and dopamine,
converting ingested? tryptophan into xx% of the brain’s serotonin.
Moreover, short-chain fatty acids (SCFAs) produced by the gut microbiome
may impact the brain directly by modulating neurotrophic factors and neuroinflammation
(Dalile et al. 2019).
Other indirect pathways for gut microbial influence on the brain include vagus nerve stimulation,
neuroendocrine modulation and immune system regulation
(Rhee, Pothoulakis, and Mayer 2009).


Rapidly growing literature connects the metabolic potential of the gut microbiome and brain function in humans
(CITE),
but the overwhelming majority of this research is performed in adults.
Importantly, both the gut microbiome and the brain undergo dramatic and rapid development over the first postnatal years
(Bonham et al. 2023; Yassour et al. 2016).
However, very little is currently known about how gut-brain influences emerge
or change during this critical window for both systems
(CITE).
Interrogating this early co-development in humans is therefore key to both understanding adaptive gut-brain function
and behavior and informing strategies to support it.
The visual cortex has been shown to be sensitive to gut microbiome modulations in adults
(CITES),
yet the visual cortex undergoes its most rapid period of plasticity and maturation over infancy
when the microbiome changes most significantly
(CITES).
The changes can be robustly indexed via EEG with the visual-evoked potential (VEP) response to visual stimuli from birth.
The VEP is especially useful for indexing visual neurodevelopment as its morphology includes amplitude deflections
and latencies to those deflections that reflect maturation of function and structure, respectively.

Here, we investigated the co-development of microbial gene functional potential
- specifically genes encoding enzymes that metabolize neuroactive compounds -
and visual neurodevelopment as indexed by the VEP in a longitudinal community samples
of 185 infants from Gugulethu in Cape Town, South Africa recruited
as part of an ongoing prospective study, called “Khula”
(Zieff et al. 2024).
Stool samples and EEG were collected at up to 3 visits in the first 18 months of life
(Figure 1A, B, Table 1;
visit 1, N = 94, age 3.7 ± 0.77 months,
visit 2, N = 77, age 8.64 ± 1.39 months,
visit 3, N = 69, age 14.2 ± 0.99 months).
Shotgun metagenomic sequencing was used to obtain microbial gene sequences from infant stool samples.
Latencies and amplitudes were extracted from each deflection component of the VEP
(i.e., first negative-going deflection, N1; first positive-going deflection, P1; and second negative-going deflection, N2),
producing six VEP features of interest.
We evaluated the concurrent association between microbial genes and VEP amplitudes and latencies
and tested prospective influences of microbial genes from early visits on VEP changes at later visits
to reveal the temporal dynamics of gut-brain co-development within individuals during this most critical window of maturation.

= Results

Here are some results.

#figure(caption: [Figure 1 - A caption])[#image(width: 50%, "placeholder.png")]

And some more...

= Discussion
#lorem(50)

= Methods
#lorem(40)

#bibliography("refs.bib", style: "nature")
