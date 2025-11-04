## IsoformCheck

Protein isoform analysis from de novo genome assemblies.

### Overview

IsoformCheck takes in haplotype phased assemblies and calls protein isoforms based on the sample sequence.
Reference transcript annotations are lifted over to samples with [liftoff](https://github.com/agshumate/Liftoff).
The coding sequences of the lifted over annotations are used to call amino acid sequences of the transcripts.
The isoforms can then be compared between samples to find novel sequences and copy count variation.

### Glossary

![Figure of sample, allele, isoform](docs/sample_allele_isoform.png)

- Sample: an assembly consisting of two haplotypes of a single individual.
- Transcript: a reference transcript, ie. a splice variant of a gene. Only protein coding genes are used.
- Isoform: a distinct amino acid sequence of a transcript.
- Allele: a single instance of an isoform in a sample.
- Allele set: all alleles of a single sample for a single transcript.
- Group: an arbitrary group of samples used for creating contingency tables and running chi squared tests.

### Detecting novel isoforms and allele sets

![Figure of comparing novel samples to a database](docs/compare_sapmles_to_database.png)

IsoformCheck can be used to detect novel isoforms and allele sets in samples.
This requires haplotype phased assemblies of the novel samples, and compares the protein coding transcripts present in the novel samples with existing samples in an IsoformCheck database.
To do this, you first need a pre-existing database of haplotype phased samples (see #prebuilt-isoformcheck-databases-of-long-read-samples).
Then use the `comparesamples` command of IsoformCheck.

The output will have three tables:
- `output_allelesets.tsv`: All allele sets of the given novel samples.
- `output_novel_isoforms.tsv`: Isoforms which are novel to the samples and not previously present in the database.
- `output_novel_allelesets.tsv`: Allele sets which are novel to the samples and not previously present in the database.

### Prebuilt IsoformCheck databases of long read samples

A database of 231 samples from the [Human Pangenome Reference Consortium](https://humanpangenome.org/) and 61 samples from the [Human Genome Structural Variation Consortium](https://www.hgsvc.org/).

### Commands

- `initialize`: Create a new IsoformCheck database. Requires a reference genome and reference annotation.
- `liftover`: Lift over annotations to a single haplotype. Can be used to parallelize the sample adding step.
- `addsample`: Adds a single haplotype to an IsoformCheck database.
- `addsamples`: Adds multiple samples to an IsoformCheck database.
- `listsamples`: List all samples in the database.
- `addgroup`: Add a sample into a group.
- `removegroup`: Remove a sample from a group.
- `listgroups:` List all groups.
- `rename`: Renames isoforms based on their coverages. Must be called after all samples have been added and before doing any analysis.
- `comparesamples`: Compares new samples to existing samples.
- `exportisoforms`: Export isoform names, sequences, and copy counts.
- `exportallelesets`: Export allele sets.
- `contingencytable`: Export contingency tables of allele sets by group.
- `chisquare`: Run chi squares tests to check which transcripts vary between groups.
- `stats`: Basic statistics about the database.
- `validate`: Check database validity.
