# MetaAHCB
MetaAHCB, short for Metagenomic Assembly, Hybrid Correction, and Binning, is a bioinformatics pipeline designed for metagenomic DNA sequencing data. It integrates genome assembly, hybrid error correction using both short- and long-read sequencing data, and contig binning into a streamlined workflow.

**Metagenomic Assembly, Hybrid Correction, and Binning**



## Project Information

**Institution:** Qingdao Institute of BGI Genomics  
**Main developer:** Jun Xia  
**Contact:** xiajun.kyu@gmail.com  

**Conceptual guidance and technical support:**  
Jianwei Chen, Mengyang Xu, Ying Sun, and Yanwei Qi  

**Version:** 1.2




## Overview

With the rapid development of third-generation sequencing (TGS) technologies, long-read sequencing has greatly improved genome assembly by providing extended read lengths. However, its relatively high sequencing error rate remains a major challenge. In contrast, second-generation sequencing (NGS) provides highly accurate short reads but is limited by short read lengths, which often results in fragmented assemblies. Combining the complementary advantages of short- and long-read sequencing is therefore essential for generating accurate and complete metagenomic assemblies.

MetaAHCB integrates hybrid error correction, metagenomic assembly, and contig binning into an automated workflow. By utilizing short reads to improve the accuracy of long reads and leveraging long reads to enhance assembly continuity, MetaAHCB aims to generate longer and more reliable contigs and high-quality MAGs from complex microbial communities.

This pipeline was developed for environmental metagenomic studies, where microbial communities contain diverse and uncultivated microorganisms. By improving genome reconstruction from metagenomic DNA sequences, MetaAHCB facilitates downstream analyses of microbial diversity, evolution, functional potential, and ecological roles.




## Program Files

After copying the MetaAHCB directory to the working directory (`./`), the core program files are organized as follows:

| File | Description |
|---|---|
| `./MetaAHCB/pipeline.sh` | Main pipeline script written in Bash shell |
| `./MetaAHCB/scripts/assign_reads_clusters.py` | Python script for assigning TGS reads into clusters |
| `./MetaAHCB/scripts/bgm_cluster.py` | Python script for machine-learning-based clustering of TGS reads |


## Workflow

The overall workflow of MetaAHCB is shown below:  

![MetaAHCB workflow](MetaAHCB_workflow.pdf)


## Pipeline Steps

**(1) Long-read preprocessing and initial assembly**

Long reads (>1000 bp) from third-generation sequencing are extracted using length filtering. The k-mer frequency profiles of individual reads are calculated, followed by machine-learning-based clustering and dimensionality reduction. Reads within each cluster are assembled independently and subsequently merged and dereplicated to generate initial contigs (contigs①).

**(2) Hybrid error correction**

Quality-controlled second-generation paired-end short reads and original long reads are sequentially used to correct contigs①. Short reads improve nucleotide-level accuracy, while long reads further correct assembly errors, generating corrected contigs (contigs②).

**(3) Short-read assembly and contig integration**

Short reads are independently assembled to generate short contigs (contigs③). These short contigs are merged with contigs② to improve assembly completeness, producing integrated contigs (contigs④). Short contigs that cannot be merged with long-read-based assemblies are retained separately as contigs⑤.

**(4) Coverage estimation**

Short reads from each sample are mapped independently to contigs④ and contigs⑤. Mapping depth information is calculated for downstream genome binning.

**(5) Metagenomic binning**

Contigs④ and contigs⑤ are separately binned using coverage information across samples, generating two sets of genome bins (bins① and bins②).

**(6) Bin dereplication**

Bins① and bins② are combined and dereplicated based on 99% average nucleotide identity (ANI). Shorter redundant genomes are removed, resulting in the final candidate bins (bins③).

**(7) MAG quality assessment and taxonomy assignment**

Candidate bins are evaluated using CheckM2 for completeness and contamination estimation and classified using GTDB-Tk. The final outputs include MAGs with different quality levels and discarded bins that do not meet quality requirements.



## Functions and Software Dependencies

MetaAHCB integrates the following tools and scripts:

| Function | Tool |
|---|---|
| Read length filtering | seqkit |
| Read quality control | fastp |
| k-mer frequency profiling | gc_nmer |
| Read clustering based on k-mer profiles | bgm_cluster.py, assign_reads_clusters.py |
| Long-read assembly | metaFlye |
| Contig merging | QuickMerge |
| Sequence dereplication | MMseqs2 |
| Short-read-based correction | FMLRC2 (HERO) |
| Long-read-based correction | HERO |
| Short-read assembly | MEGAHIT |
| Read mapping | minimap2 |
| SAM/BAM processing | samtools |
| Coverage calculation | jgi_summarize_bam_contig_depths |
| Multi-sample depth integration | merge_depths.pl |
| Genome binning | MetaBAT2 |
| Contig classification | seqkit |
| Bin dereplication | dRep |
| Genome quality assessment | CheckM2 |
| Taxonomic classification | GTDB-Tk |




## Performance

The computational performance of MetaAHCB depends on the sequencing depth, number of samples, and available computing resources. Representative benchmark tests are shown below.

**Example 1: Mock community dataset**

A test dataset containing five prokaryotic species was analyzed, with approximately 5×–30× sequencing depth for each species. Using 8 CPU cores, the complete workflow required approximately 24 hours, with a peak memory usage of ~160 GB.

**Example 2: Single metagenomic sample**

A single metagenomic sample containing 10 GB of third-generation sequencing data and 20 GB of second-generation sequencing data was processed. Using 16 CPU cores, the complete workflow required approximately 96 hours, with a peak memory usage of ~500 GB. When the long-read correction step using HERO was skipped, the runtime was reduced to approximately 48 hours, with a peak memory usage of ~300 GB.

**Example 3: Multiple metagenomic samples**

Two groups, each containing five metagenomic samples, were processed. Each sample included 10 GB of third-generation sequencing data and 20 GB of second-generation sequencing data. When running `pipeline.sh` for each group with 16 CPU cores, the complete workflow required approximately 20 days, with a peak memory usage of ~1000 GB. Without the HERO long-read correction step, the runtime was reduced to approximately 10 days, with a peak memory usage of ~600 GB.

(Performance optimization and customization of individual workflow steps are not described in detail here.)




## Software Requirements

MetaAHCB requires the following software and libraries:

### Python environment

- Python 3
- Python packages:
  - gzip
  - pandas
  - Biopython (`Bio`)
  - argparse
  - threading
  - os
  - numpy
  - scikit-learn (`sklearn`)



## External software

- fastp
- seqkit
- gc_nmer
- metaFlye
- MMseqs2
- minimap2
- HERO
- samtools
- MEGAHIT
- QuickMerge
- MetaBAT2
- dRep
- GTDB-Tk
- CheckM2
