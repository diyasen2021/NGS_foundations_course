# Module 1: Introduction to Genomics and Sequencing Technologies

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 1 | Module 1 of 4

---

## Table of Contents

- [1.1 What is Genomics?](#11-what-is-genomics)
- [1.2 The Central Dogma of Molecular Biology](#12-the-central-dogma-of-molecular-biology)
- [1.3 Structure of the Human Genome](#13-structure-of-the-human-genome)
- [1.4 From Sanger Sequencing to Next Generation Sequencing](#14-from-sanger-sequencing-to-next-generation-sequencing)
- [1.5 Major NGS Platforms](#15-major-ngs-platforms)
- [1.6 Choosing the Right Sequencing Technology](#16-choosing-the-right-sequencing-technology)
- [1.7 Key Terminology](#17-key-terminology)

---

## 1.1 What is Genomics?

Genomics is the branch of molecular biology that studies the complete set of DNA — the genome — of an organism. Unlike classical genetics, which focuses on individual genes and their roles in inheritance, genomics takes a broader view, examining the structure, function, evolution, and mapping of all genes and non-coding sequences together. The human genome contains approximately **3.2 billion base pairs** of DNA packed into 23 pairs of chromosomes, encoding roughly **20,000–25,000 protein-coding genes**. Remarkably, these genes account for only about 2% of the total genome; the remaining 98% consists of regulatory elements, introns, repetitive sequences, and regions whose functions are still being actively investigated.

Genomics has transformed biological research by shifting the focus from studying one gene at a time to interrogating entire genomes simultaneously. This shift has been made possible largely by advances in DNA sequencing technologies, which have become faster, cheaper, and more accurate over the past four decades.

---

## 1.2 The Central Dogma of Molecular Biology

To understand why sequencing matters, students must first be comfortable with the Central Dogma, proposed by Francis Crick in 1958. It describes the flow of genetic information within a biological system:

$$\text{DNA} \rightarrow \text{RNA} \rightarrow \text{Protein}$$

DNA is transcribed into messenger RNA (mRNA) by the enzyme RNA polymerase. The mRNA is then translated into a protein by ribosomes, with transfer RNAs (tRNAs) reading the codons in triplets. Each codon corresponds to a specific amino acid, and the sequence of amino acids ultimately determines the structure and function of the protein.

A few important nuances for students to appreciate:

- **Not all DNA is transcribed, and not all RNA is translated.** Regulatory RNAs such as miRNA, lncRNA, and rRNA play important structural and regulatory roles without ever becoming protein.
- **Reverse transcription** — RNA being converted back to DNA — occurs in retroviruses (like HIV) and is also harnessed in the laboratory for RNA sequencing (RNA-seq), where mRNA is converted to complementary DNA (cDNA) before sequencing.
- **Epigenetic modifications** (methylation of DNA, histone modification) can alter gene expression without changing the underlying DNA sequence — an important concept when interpreting NGS data from epigenomics experiments.

---

## 1.3 Structure of the Human Genome

The genome is organised hierarchically. At the finest level, DNA is a double-stranded helix made up of four nucleotide bases — **Adenine (A)**, **Thymine (T)**, **Guanine (G)**, and **Cytosine (C)** — paired according to complementary base-pairing rules (A–T and G–C). These strands are read in a **5' to 3' direction**, and this directionality is critical for understanding how sequencing reads are generated and mapped.

At a higher level, the genome is divided into:

**Exons**  
The coding portions of a gene that are retained in the mature mRNA after splicing and ultimately translated into protein. For a typical human gene, exons make up a small fraction of the total gene length.

**Introns**  
Non-coding sequences that interrupt the exons within a gene. They are transcribed into pre-mRNA but are spliced out before translation. Introns were long considered "junk DNA," but research has revealed that many harbour regulatory sequences and small non-coding RNAs.

**Promoters and Regulatory Regions**  
Sequences upstream (and sometimes downstream) of a gene that control when, where, and how strongly a gene is expressed. Transcription factors bind these regions to activate or repress transcription.

**Repetitive Elements**  
Sequences such as SINEs (Short Interspersed Nuclear Elements), LINEs (Long Interspersed Nuclear Elements), and satellite repeats make up over 50% of the human genome. They pose a significant challenge in NGS analysis because short reads from these regions are difficult to map uniquely.

**Non-coding RNAs**  
Includes ribosomal RNA (rRNA), transfer RNA (tRNA), microRNA (miRNA), and long non-coding RNA (lncRNA) — all transcribed but not translated.

> 💡 **Practical note:** Understanding genome structure is directly relevant to experiment design. If you are only interested in protein-coding regions, you might choose **Whole Exome Sequencing (WES)** rather than **Whole Genome Sequencing (WGS)**, reducing both cost and data complexity considerably.

---

## 1.4 From Sanger Sequencing to Next Generation Sequencing

### Sanger Sequencing (First Generation)

The story of DNA sequencing begins with Frederick Sanger, who developed chain-termination sequencing in **1977** — a method that earned him his second Nobel Prize. The principle relies on incorporating **dideoxynucleotides (ddNTPs)** — modified bases that lack the 3'-OH group required for chain elongation. 
- In Sanger sequencing, special bases (Dideoxynucleotides) stop DNA synthesis at random positions, producing fragments of different lengths.
- Each fragment tells you the identity of the last base
- By ordering fragments by size (through gel electrophoresis), the DNA sequence can be reconstructed 

Modern automated Sanger sequencing replaced the radioactive gels with fluorescently labelled ddNTPs and capillary electrophoresis, allowing sequencing of a single fragment up to about **600–1000 base pairs** in a single run.

Sanger sequencing is still used today for:

- Validating variants discovered by NGS
- Sequencing short, targeted regions (e.g. checking a PCR product)
- Clinical confirmation of pathogenic mutations before reporting

Its fundamental limitation is **throughput** — it can only sequence one fragment at a time, making it impractical for sequencing entire genomes at scale.

### The NGS Revolution

Next Generation Sequencing (also called second-generation or massively parallel sequencing) emerged in the mid-2000s and fundamentally changed the field. Rather than sequencing one DNA fragment at a time, NGS sequences **millions of fragments simultaneously** in a single instrument run. This massive parallelisation reduced the cost of sequencing a human genome from approximately **$3 billion** (Human Genome Project, completed 2003) to under **$1,000** today — a cost reduction faster than Moore's Law.

- The key conceptual difference is that NGS breaks the DNA into many small fragments, sequences them all in parallel, and then uses computational alignment to reassemble the sequence by mapping each short read back to a reference genome.
- This introduces two important trade-offs: reads are much shorter than Sanger reads (typically **75–300 bp** for Illumina), and the analysis requires significant bioinformatics expertise.

---

## 1.5 Major NGS Platforms

### Illumina — Sequencing by Synthesis (SBS)

Illumina is by far the most widely used NGS platform globally, dominating both research and clinical sequencing. 
- Its chemistry is based on **reversible terminator sequencing by synthesis**.
- DNA fragments are attached to a flow cell surface and amplified locally to form dense clusters.
- During each sequencing cycle, fluorescently labelled reversible terminator nucleotides are incorporated one at a time.
- The fluorescence is imaged, the terminator is chemically removed, and the cycle repeats.
- This generates highly accurate short reads (75–300 bp) with error rates as low as **0.1%**.

Illumina instruments range from the small **MiSeq** (suitable for small genomes or targeted panels) to the ultra-high-throughput **NovaSeq X Plus**, capable of sequencing over 20,000 whole human genomes per year. For most UG and PG projects involving WGS, WES, RNA-seq, or ChIP-seq, Illumina data will be the primary data type encountered.

### Oxford Nanopore Technologies (ONT)

Nanopore sequencing represents a fundamentally different approach. DNA or RNA is passed through a biological nanopore (a protein channel) embedded in a membrane. As individual bases move through the pore, they cause characteristic disruptions to an electrical current flowing across the membrane. These current signatures are decoded by a neural-network-based base caller into nucleotide sequences.

The major advantage of nanopore sequencing is **read length** — reads can span from a few kilobases to over **4 megabases**, making it exceptional for resolving repetitive regions, structural variants, and complete chromosome-scale assemblies. The **MinION** device is roughly the size of a USB stick and has been used in field settings and even on the International Space Station. The main limitation has historically been a higher raw error rate (~5–10%), though recent chemistry and basecalling improvements have substantially closed this gap.

### PacBio — Single Molecule Real Time (SMRT) Sequencing

PacBio uses a **zero-mode waveguide** approach, where individual DNA polymerase molecules sit at the bottom of tiny wells and synthesise DNA in real time. Each nucleotide is fluorescently labelled, and the signal is detected as it is incorporated. Like Nanopore, PacBio generates long reads — typically **15–25 kb** for standard long-read mode.

The newer **HiFi (CCS) mode** sequences the same circular DNA molecule multiple times, generating consensus reads with accuracy exceeding **99.9%**, combining the length of long reads with the accuracy approaching short reads. PacBio HiFi is increasingly used for phased genome assembly, structural variant detection, and full-length isoform sequencing (Iso-Seq).

### Ion Torrent

Ion Torrent sequencing detects **hydrogen ions** released during DNA synthesis rather than fluorescence. As a nucleotide is incorporated, a proton is released, causing a measurable change in pH that is detected by a semiconductor chip. The technology is faster and less expensive to instrument than optical methods but generates shorter reads (200–600 bp) with a higher error rate in homopolymer regions. Ion Torrent is commonly used in clinical settings for targeted sequencing panels, such as oncology panels testing a defined set of cancer-associated variants.

### Platform Comparison Summary

| Platform | Read Length | Accuracy | Throughput | Best For |
|---|---|---|---|---|
| Illumina | 75–300 bp | ~99.9% | Up to 6 Tb | WGS, RNA-seq, ChIP-seq |
| Oxford Nanopore | 1 kb – 4 Mb | ~90–99% | Up to 50 Gb | Structural variants, metagenomics |
| PacBio HiFi | 15–25 kb | ~99.9% | ~360 Gb | Genome assembly, phasing |
| Ion Torrent | 200–600 bp | ~99% | ~15 Gb | Targeted panels, amplicon |

---

## 1.6 Choosing the Right Sequencing Technology

One of the most important practical skills for any genomics researcher is matching the sequencing technology to the biological question. There is no single "best" platform — the right choice depends on read length requirements, accuracy needs, budget, sample type, and downstream analysis.

| Research Question | Recommended Platform |
|---|---|
| Routine variant calling (SNPs, small indels) | Illumina — best accuracy, cost, and coverage |
| De novo genome assembly | PacBio HiFi or ONT long reads |
| Structural variant detection | ONT or PacBio |
| Portable / field sequencing | Oxford Nanopore MinION |
| Clinical targeted panels | Ion Torrent or Illumina MiSeq |
| Full-length isoform sequencing | PacBio Iso-Seq or ONT cDNA |

> ⚠️ **Important:** Choosing the wrong platform at the start of a project can result in data that cannot answer your biological question. Always define your question first, then select the technology.

---

## 1.7 Key Terminology

Before moving into the practical modules, students should be familiar with the following terms:

| Term | Definition |
|---|---|
| **Read** | A single sequenced DNA fragment, represented as a string of A, T, G, C characters |
| **Read Length** | The number of base pairs in a single read (e.g. Illumina: 150 bp paired-end) |
| **Paired-end Sequencing** | Sequencing both ends of a DNA fragment — improves alignment and variant detection |
| **Coverage (Depth)** | Average number of reads covering each genome position. 30× = each base covered ~30 times |
| **FASTQ Format** | Standard raw data file format — contains sequence + per-base quality scores |
| **Phred Quality Score** | Log-scale accuracy measure. Q30 = 99.9% accuracy (1 in 1,000 error chance) |
| **Reference Genome** | Representative DNA sequence used as a standard for alignment (human: GRCh38/hg38) |
| **Genome Assembly** | Reconstructing genome sequence from reads — reference-guided or de novo |
| **WGS** | Whole Genome Sequencing — sequences the entire genome |
| **WES** | Whole Exome Sequencing — sequences only protein-coding regions (~2% of genome) |

---

> *This foundation in genomics and sequencing technologies sets the stage for everything that follows in the workshop. A solid grasp of these concepts — the structure of the genome, the principles behind each sequencing platform, and the vocabulary of NGS data — will make the hands-on bioinformatics sessions considerably more intuitive.*

---
