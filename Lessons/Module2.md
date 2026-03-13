# Module 2: NGS Workflow and Data Generation

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 1 | Module 2 of 4

---

## Table of Contents

- [2.1 Overview of the NGS Workflow](#21-overview-of-the-ngs-workflow)
- [2.2 Step 1 — Sample Collection and Nucleic Acid Extraction](#22-step-1--sample-collection-and-nucleic-acid-extraction)
- [2.3 Step 2 — Library Preparation](#23-step-2--library-preparation)
- [2.4 Step 3 — Cluster Generation or Loading](#24-step-3--cluster-generation-or-loading)
- [2.5 Step 4 — Sequencing](#25-step-4--sequencing)
- [2.6 Step 5 — Base Calling and Quality Scoring](#26-step-5--base-calling-and-quality-scoring)
- [2.7 The FASTQ File Format](#27-the-fastq-file-format)
- [2.8 Types of NGS Experiments](#28-types-of-ngs-experiments)
- [2.9 Sequencing Depth and Coverage — How Much Data Do You Need?](#29-sequencing-depth-and-coverage--how-much-data-do-you-need)
- [2.10 Common Sources of Error and Bias](#210-common-sources-of-error-and-bias)

---

## 2.1 Overview of the NGS Workflow

Regardless of the sequencing platform or biological question, all NGS experiments follow the same fundamental series of steps. Understanding this end-to-end workflow is essential before touching any data, because decisions made at each stage have downstream consequences that cannot be undone computationally.

The core workflow is:

```
Sample Collection & QC
        ↓
Nucleic Acid Extraction
        ↓
Library Preparation
        ↓
Cluster Generation / Flow Cell Loading
        ↓
Sequencing
        ↓
Base Calling
        ↓
FASTQ Files  ←── Your starting point for bioinformatics
```

Each of these steps introduces its own opportunities for error, bias, and variation. A poor-quality sample going into the library preparation will produce poor-quality data coming out of the sequencer — no amount of computational processing can recover biological signal that was never captured in the first place.

---

## 2.2 Step 1 — Sample Collection and Nucleic Acid Extraction

### Sample Types

The starting material for an NGS experiment can be almost any biological specimen, but the most common sources encountered in research and clinical settings include:

- **Whole blood** — rich in white blood cells; easy to collect and process; the workhorse of germline genomics
- **Tissue biopsies** — fresh-frozen tissue preserves nucleic acids well; FFPE (formalin-fixed paraffin-embedded) tissue is common in clinical oncology but introduces DNA damage
- **Saliva and buccal swabs** — non-invasive; widely used in population studies
- **Cell lines** — well-characterised; useful for benchmarking and functional studies
- **Environmental samples** — soil, water, gut content; used in metagenomics
- **FFPE tissue** — common in clinical archives but introduces formalin-induced cross-links and DNA fragmentation that complicate sequencing

### DNA Extraction

For DNA sequencing (WGS, WES, targeted panels, ChIP-seq), the goal is to extract high-molecular-weight, high-purity genomic DNA. Common methods include:

- **Column-based kits** (e.g. Qiagen DNeasy) — fast and convenient; suitable for most sample types
- **Phenol-chloroform extraction** — labour-intensive but yields high-quality, high-molecular-weight DNA; important for long-read sequencing where fragment length matters
- **Magnetic bead-based methods** — high-throughput; widely used in automated pipelines

### RNA Extraction

For RNA sequencing (RNA-seq), extraction must preserve RNA integrity. RNA is inherently less stable than DNA due to the ubiquitous presence of RNases.

- **TRIzol / TRI Reagent** — a monophasic solution that simultaneously disrupts cells and separates RNA from DNA and proteins
- **Column-based kits** (e.g. Qiagen RNeasy) — convenient and reproducible
- **RNA later** — a stabilisation reagent added immediately after collection to protect RNA from degradation during storage and transport

### Quality Control of Extracted Nucleic Acids

Before proceeding to library preparation, nucleic acid quality and quantity must be assessed. Using poor-quality input material is one of the most common and costly mistakes in an NGS experiment.

| QC Metric | Tool | What to Check |
|---|---|---|
| Concentration | Qubit fluorometer | Accurate quantification (spectrophotometers overestimate due to free nucleotides) |
| Purity | NanoDrop (A260/A280, A260/A230) | A260/A280 ~1.8 for DNA, ~2.0 for RNA; lower values indicate protein or solvent contamination |
| Integrity (DNA) | Agarose gel or TapeStation | High-MW band with no smearing; important for long-read sequencing |
| Integrity (RNA) | Bioanalyzer or TapeStation | RIN (RNA Integrity Number) ≥ 7 recommended; ≥ 8 for gene expression studies |

> 💡 **Tip for Indian research settings:** Qubit is strongly preferred over NanoDrop for concentration measurement for NGS input. NanoDrop reads at 260 nm and cannot distinguish intact nucleic acids from free nucleotides, degradation products, or contaminants — all of which absorb at the same wavelength.

---

## 2.3 Step 2 — Library Preparation

Library preparation is the process of converting your extracted nucleic acid into a form that the sequencer can read. This is arguably the most critical and variable step in the entire workflow, and it is where most experimental bias is introduced.

### What is a Sequencing Library?

A sequencing library is a collection of DNA fragments of appropriate size, each flanked by **adapter sequences** that allow the fragments to:

1. Bind to the flow cell surface (Illumina) or be loaded into a nanopore
2. Be amplified (where required)
3. Be identified by the sequencer during imaging or signal detection
4. Be demultiplexed if multiple samples are pooled together (multiplexing)

### General Steps in DNA Library Preparation

**1. Fragmentation**  
Genomic DNA must be broken into fragments of the appropriate size for the sequencing platform. For Illumina, the optimal insert size is typically **150–500 bp**. Fragmentation is achieved by:
- **Sonication** (e.g. Covaris) — physical shearing using focused ultrasound; highly reproducible
- **Enzymatic fragmentation** — uses restriction enzymes or non-specific nucleases; faster but can introduce sequence bias

**2. End Repair**  
Fragmentation leaves ragged, non-blunt ends. End repair uses a combination of DNA polymerase and exonuclease activity to generate blunt-ended, 5'-phosphorylated fragments.

**3. A-tailing**  
A single adenosine (A) base is added to the 3' end of each blunt fragment by a polymerase. This A-overhang is complementary to a T-overhang on the sequencing adapters, enabling efficient ligation.

**4. Adapter Ligation**  
Platform-specific adapter sequences are ligated to both ends of the fragment. These adapters contain:
- Sequencing primer binding sites
- Flow cell binding sequences (for Illumina)
- Index (barcode) sequences for multiplexing

**5. Size Selection**  
AMPure XP magnetic beads or gel-based size selection remove adapter dimers (adapters that have ligated to each other without an insert), very small fragments, and very large fragments outside the target range.

**6. PCR Amplification**  
Most standard library protocols include a PCR amplification step to increase the amount of adapter-ligated material. Typically 8–15 cycles are used. Key concerns:
- Too many PCR cycles introduce **PCR duplicates** — multiple copies of the same original molecule that can be mistaken for independent evidence of a variant
- PCR amplification can also introduce GC bias, preferentially amplifying AT-rich or GC-rich regions less efficiently

> 💡 **PCR-free libraries** are available for high-coverage WGS and are strongly preferred when studying difficult regions of the genome (e.g. highly repetitive regions, GC-extremes). They reduce duplicate rates and GC bias at the cost of requiring more input DNA (typically ≥ 1 µg).

### RNA-seq Library Preparation

For RNA-seq, there is an additional step before library preparation: the RNA must either be **enriched for mRNA** or **depleted of ribosomal RNA (rRNA)**, which otherwise accounts for 80–90% of total cellular RNA.

- **Poly-A selection** — uses oligo-dT beads to capture polyadenylated mRNA; simple and cost-effective; misses non-polyadenylated transcripts
- **Ribosomal RNA depletion** — removes rRNA using complementary probes; preserves lncRNA, pre-mRNA, and other non-polyadenylated transcripts; better for degraded samples where poly-A tails may be lost

After enrichment or depletion, the RNA is fragmented, converted to cDNA by reverse transcriptase, and then library preparation proceeds similarly to DNA.

**Strand-specific (stranded) libraries** preserve information about which DNA strand was originally transcribed, which is critical for accurately quantifying overlapping genes on opposite strands. Most modern RNA-seq protocols use stranded library preparation.

### Multiplexing and Index Sequences

Running a single sample per sequencing lane is expensive and wasteful. In practice, multiple samples are pooled together in a single lane by adding a short unique **index sequence** (also called a barcode) to the adapters of each sample's library during preparation. After sequencing, the reads are **demultiplexed** — computationally separated back into individual sample files — based on their index sequence.

Dual indexing (i7 + i5) uses two independent barcodes per sample and is now standard practice, as it dramatically reduces index hopping artefacts that can occur on patterned flow cells.

---

## 2.4 Step 3 — Cluster Generation or Loading

### Illumina — Cluster Generation

On Illumina platforms, the library is loaded onto a **flow cell** — a glass slide with millions of tiny wells (patterned flow cells, e.g. NovaSeq) or a lawn of oligonucleotides (older random flow cells). The adapter sequences on the library fragments hybridise to complementary oligos on the flow cell surface.

Each bound fragment is then amplified locally through a process called **bridge amplification**, in which the fragment bends over and hybridises to a nearby surface oligo, and is extended. This is repeated many times, generating a **cluster** of ~1,000 identical copies of the original molecule. Clusters are necessary because the fluorescent signal from a single molecule would be too weak to detect reliably.

On patterned flow cells (HiSeq X, NovaSeq), fragments are loaded into pre-formed nanowells at a density high enough that most wells contain exactly one cluster, reducing the problem of mixed clusters (polyclonal clusters that reduce data quality).

### Oxford Nanopore — Flow Cell Loading

For nanopore sequencing, there is no amplification. Native DNA or RNA molecules are loaded directly onto a flow cell containing an array of nanopores embedded in a membrane. A motor protein controls the speed at which each molecule threads through the pore. Because nanopore sequencing is amplification-free, it can directly detect base modifications such as 5-methylcytosine without bisulfite conversion — a significant advantage for epigenomics.

---

## 2.5 Step 4 — Sequencing

### Illumina Sequencing by Synthesis

Once clusters are generated, sequencing begins. In each cycle:

1. All four fluorescently labelled, reversible-terminator nucleotides are washed over the flow cell
2. Each cluster incorporates exactly one nucleotide (determined by the next base in the template)
3. Unincorporated nucleotides are washed away
4. The flow cell is imaged — the colour of fluorescence identifies the base
5. The terminator group and fluorescent dye are chemically cleaved, regenerating a free 3'-OH
6. The cycle repeats

This process continues for the full read length — typically **150 cycles per read** for paired-end 150 bp sequencing. In paired-end sequencing, the flow cell is flipped after completing the first read, and sequencing proceeds from the other end of the fragment for a second read of equal length.

### Sequencing Run Output

A typical Illumina run produces raw data in **BCL (Base Call) format**, which stores the intensity and base call for each cycle, for each cluster, on each tile of the flow cell. BCL files are then converted to **FASTQ format** using the `bcl2fastq` or `BCLConvert` software, and samples are demultiplexed by their index sequences simultaneously.

---

## 2.6 Step 5 — Base Calling and Quality Scoring

### What is Base Calling?

Base calling is the process of converting the raw signal from the sequencer into a nucleotide sequence. For Illumina, this means interpreting the fluorescent intensity images from each cycle. For Oxford Nanopore, it means interpreting the raw electrical current trace (called a squiggle) for each molecule passing through a pore — a considerably more complex signal processing problem solved by deep learning models (e.g. Guppy, Dorado).

### Phred Quality Scores

Every base in an NGS read is assigned a **Phred quality score (Q score)** that represents the base caller's confidence in that base call. The Q score is defined as:

$$Q = -10 \cdot \log_{10}(P)$$

Where *P* is the estimated probability that the base call is incorrect. This gives us:

| Phred Score | Error Probability | Accuracy |
|---|---|---|
| Q10 | 1 in 10 | 90% |
| Q20 | 1 in 100 | 99% |
| Q30 | 1 in 1,000 | 99.9% |
| Q40 | 1 in 10,000 | 99.99% |

A commonly used threshold is **Q30**, meaning that at least 75–80% of bases in a run should have Q ≥ 30. This is often reported as a run quality metric and is one of the first things to check when evaluating sequencing data.

### Why Quality Scores Decrease Towards the End of Reads

On Illumina platforms, base quality tends to decline towards the 3' end of reads. This is because:

- Over many cycles, clusters become slightly out-of-phase (some molecules in a cluster are one cycle ahead or behind), causing signal blurring
- Fluorescent dyes and reversible terminators degrade over the course of the run
- Enzyme efficiency decreases over time

This is why **read trimming** (covered in Module 4) is an important preprocessing step — removing low-quality bases from the ends of reads before alignment significantly improves downstream analysis accuracy.

---

## 2.7 The FASTQ File Format

The FASTQ file is the universal starting point for all NGS bioinformatics analysis. Every read produced by the sequencer is represented as four lines in the FASTQ file:

```
@SEQ_ID
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
+
BBBFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFF
```

| Line | Content | Description |
|---|---|---|
| Line 1 | `@SEQ_ID` | Read identifier — begins with `@`; contains instrument, run, flowcell, lane, tile, cluster coordinates, and pair information |
| Line 2 | Sequence | The nucleotide sequence of the read (A, T, G, C, or N for undetermined) |
| Line 3 | `+` | Separator line — sometimes repeats the read ID, often just `+` |
| Line 4 | Quality string | Per-base quality scores encoded as ASCII characters |

### Decoding the Quality String

Quality scores are encoded as ASCII characters to keep the file compact. The encoding used in modern data (Phred+33, also called Sanger encoding) maps Q scores to ASCII as:

```
Q score + 33 = ASCII character code
```

So Q30 → ASCII 63 → `?`, Q40 → ASCII 73 → `I`, and so on. The lowest printable character `!` represents Q0, while `~` represents Q93.

### FASTQ File Sizes

FASTQ files can be very large. A single whole human genome sequenced at 30× coverage on Illumina (paired-end 150 bp) generates approximately:
- ~600 million read pairs
- ~180 Gb of uncompressed FASTQ data
- ~40–60 Gb when gzip-compressed (`.fastq.gz`)

For this reason, FASTQ files are almost always stored and transferred in gzip-compressed format, and all major bioinformatics tools accept `.fastq.gz` input directly without requiring decompression.

---

## 2.8 Types of NGS Experiments

Different biological questions require different NGS experimental designs. Below are the most commonly encountered in research:

### Whole Genome Sequencing (WGS)

Sequences the entire genome, including coding and non-coding regions. Used for:
- Germline variant discovery (SNPs, indels, structural variants, CNVs)
- De novo genome assembly
- Population genomics
- Somatic mutation profiling in cancer

Typical coverage: **30× for germline; 60–100× for somatic/cancer**

### Whole Exome Sequencing (WES)

Captures and sequences only the protein-coding exons of the genome using a targeted enrichment step (hybrid capture). Covers ~2% of the genome but ~85% of known disease-causing variants. Used for:
- Rare disease diagnosis
- Clinical genetics
- Cost-effective variant discovery when only coding regions are of interest

Typical coverage: **100× mean; ≥20× for ≥95% of target bases**

### RNA Sequencing (RNA-seq)

Sequences the transcriptome — the complete set of RNA transcripts present in a sample at a given time. Used for:
- Differential gene expression analysis
- Alternative splicing analysis
- Fusion gene detection
- Novel transcript discovery
- Single-cell transcriptomics (scRNA-seq)

Typical depth: **20–50 million read pairs per sample for bulk RNA-seq**

### ChIP-seq

Combines chromatin immunoprecipitation with sequencing to map the genome-wide binding sites of transcription factors or histone modifications. Used to study gene regulation and chromatin structure.

### ATAC-seq

Assay for Transposase-Accessible Chromatin using sequencing. Maps open chromatin regions — areas of the genome accessible to transcription factors — genome-wide. A powerful tool for studying gene regulatory landscapes.

### Targeted / Amplicon Sequencing

Sequences a predefined set of genomic regions at very high depth. Uses either hybrid capture or PCR amplification to enrich target regions. Common in:
- Clinical oncology panels (e.g. BRCA1/2, TP53, EGFR)
- Pharmacogenomics testing
- Pathogen typing and surveillance

### Bisulfite Sequencing (BS-seq / WGBS)

Treats DNA with sodium bisulfite, which converts unmethylated cytosines to uracil (read as thymine) while leaving methylated cytosines unchanged. Used to map DNA methylation genome-wide.

---

## 2.9 Sequencing Depth and Coverage — How Much Data Do You Need?

Coverage is one of the most important experimental design decisions, and it directly determines the cost of the experiment. Coverage requirements depend on the application:

| Application | Recommended Coverage |
|---|---|
| Germline WGS (SNPs, small indels) | 30× |
| Somatic WGS (tumour/normal pair) | 60× tumour / 30× normal |
| WES | 100× mean; ≥20× for ≥95% targets |
| RNA-seq (bulk, standard DE) | 20–50M read pairs |
| RNA-seq (rare isoforms, lncRNA) | 100M+ read pairs |
| ChIP-seq | 20–40M reads (TF); 40–80M (histone marks) |
| ATAC-seq | 50–200M read pairs |
| Targeted panel (clinical) | 500× – 2,000× |

### The Relationship Between Coverage, Sensitivity, and Variant Allele Frequency

Coverage is not just about statistical confidence — it fundamentally determines what variants you can detect. For somatic variant calling (e.g. detecting mutations in a tumour), a mutation present in only 5% of cells (variant allele frequency, VAF = 0.05) requires very high coverage to reliably detect. At 30× coverage, you would expect only ~1–2 reads supporting that variant — far too few to distinguish signal from sequencing noise. At 500× coverage, you would expect ~25 supporting reads, which is detectable.

This is why clinical ctDNA (liquid biopsy) panels routinely sequence at **1,000× – 10,000×** to detect rare circulating tumour DNA variants.

---

## 2.10 Common Sources of Error and Bias

Understanding where errors and biases are introduced in the NGS workflow is critical for correctly interpreting results. Many apparent biological signals in NGS data are in fact technical artefacts.

### PCR Duplicates

When the same original DNA molecule is amplified multiple times during library preparation, all copies produce identical reads after sequencing. These **PCR duplicates** are not independent observations of the genome — they represent one molecule counted multiple times. In variant calling, duplicates can artificially inflate the apparent evidence for a variant.

Tools such as **Picard MarkDuplicates** or **samtools markdup** identify and flag duplicate reads after alignment based on their identical start positions on both strands. Most variant callers ignore flagged duplicates.

### GC Bias

PCR amplification is less efficient at extreme GC content (very high or very low). This results in uneven coverage across the genome, with GC-rich or AT-rich regions having lower read depth. This is a particular problem for WGS applications requiring uniform coverage, and is one reason PCR-free libraries are preferred for high-quality WGS.

### Adapter Contamination

If a DNA insert is shorter than the read length, the sequencer will read into the adapter sequence on the other end. Adapter sequences in reads will fail to align to the reference genome and can also introduce artefacts in downstream analyses. Adapter trimming (Module 4) removes these sequences.

### Index Hopping

On Illumina patterned flow cells (HiSeq 3000/4000, NovaSeq), a small proportion of library molecules (~0.1–2%) can switch their index sequences during cluster generation. This means a small number of reads from one sample are assigned to another sample after demultiplexing. The use of **unique dual indexes (UDIs)** — where each sample has a unique combination of both i7 and i5 indexes not shared with any other sample in the pool — virtually eliminates misassignment from index hopping.

### FFPE Artefacts

Formalin fixation causes cytosine deamination, creating C→T (and G→A on the opposite strand) artefact mutations. These can mimic real somatic mutations and are a significant problem in cancer genomics using archival FFPE tissue. Specific FFPE repair kits and bioinformatic filters are used to reduce this artefact.

### Strand Bias

Certain sequencing errors occur preferentially on one strand of the DNA. If a variant is only supported by reads on one strand (forward or reverse), this is suspicious and may indicate a sequencing artefact rather than a true variant. Most variant callers report strand bias as a filter flag (e.g. `FS` — Fisher Strand in GATK).

---

> 💡 **Key takeaway:** The quality of your NGS data is determined long before you open a terminal. Every decision — sample type, extraction method, library kit, sequencing depth, and platform — shapes what you can and cannot detect. Understanding the wet-lab workflow makes you a much better computational analyst, because you can distinguish genuine biological signal from technical noise.

---

**Previous:** [Module 1 — Introduction to Genomics & Sequencing Technologies](./module1_introduction_to_genomics.md)  
**Next:** [Module 3 — Quality Control of Raw Sequencing Data](./module3_quality_control.md)
