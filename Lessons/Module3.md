# Module 3: Quality Control of Raw Sequencing Data

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 1 | Module 3 of 4 | 🧪 Hands-on

---

## Table of Contents

- [3.1 Why QC Matters](#31-why-qc-matters)
- [3.2 What Are We Looking For?](#32-what-are-we-looking-for)
- [3.3 Tool Overview](#33-tool-overview)
- [3.4 FastQC — In Depth](#34-fastqc--in-depth)
- [3.5 Interpreting FastQC Reports — Module by Module](#35-interpreting-fastqc-reports--module-by-module)
- [3.6 MultiQC — Aggregating Reports Across Samples](#36-multiqc--aggregating-reports-across-samples)
- [3.7 QC for Long-Read Data](#37-qc-for-long-read-data)
- [3.8 Hands-on Exercises](#38-hands-on-exercises)
- [3.9 QC Decision Framework — Pass, Warn, or Fail?](#39-qc-decision-framework--pass-warn-or-fail)
- [3.10 Common QC Problems and Their Causes](#310-common-qc-problems-and-their-causes)

---

## 3.1 Why QC Matters

Receiving a FASTQ file from a sequencing facility does not mean the data is ready for analysis. Raw sequencing data almost always contains:

- Bases called with low confidence, particularly towards the 3' end of reads
- Adapter sequences that were ligated during library preparation and not removed
- PCR duplicates introduced during amplification
- Reads contaminated with sequences from other organisms
- Overrepresented sequences caused by PCR amplification of specific fragments
- Systematic biases in base composition or GC content

Running quality control before any analysis is not optional — it is the first mandatory step in every NGS pipeline. Proceeding with poor-quality data leads to inaccurate alignments, false variant calls, inflated or deflated gene expression estimates, and ultimately incorrect biological conclusions.

The goal of QC at this stage is **not** to fix the data (that is the job of preprocessing, covered in Module 4) but to **understand it** — to know what you are working with, identify potential problems, and make informed decisions about how to proceed.

---

## 3.2 What Are We Looking For?

Before running any tool, it is useful to know what a high-quality dataset looks like versus a problematic one. The key metrics to evaluate are:

| Metric | What it tells you |
|---|---|
| Per-base sequence quality | Are base calls reliable across the length of the read? |
| Per-sequence quality scores | What is the distribution of average quality across all reads? |
| Per-base sequence content | Is there unexpected base composition bias at any position? |
| Per-sequence GC content | Does the GC distribution match the expected genome GC content? |
| Sequence length distribution | Are reads the expected length, or has trimming already occurred? |
| Sequence duplication levels | Is there evidence of excessive PCR amplification? |
| Adapter content | Are adapter sequences present in reads? |
| Overrepresented sequences | Are any specific sequences present far more than expected by chance? |
| N content | Are there positions where the base caller could not determine the base? |

---

## 3.3 Tool Overview

| Tool | Best For | Output |
|---|---|---|
| **FastQC** | Per-sample QC of Illumina short reads | HTML report + zip archive per FASTQ |
| **MultiQC** | Aggregating FastQC (and other tool) reports across many samples | Single HTML report |
| **FastQ Screen** | Checking for cross-species contamination | HTML report showing which genomes reads map to |
| **NanoStat** | Summary statistics for Oxford Nanopore FASTQ/FAST5 | Text/TSV summary |
| **NanoPlot** | Visualising Nanopore read length and quality distributions | HTML report with plots |
| **PycoQC** | Interactive QC for Nanopore data | Interactive HTML report |
| **Falco** | Faster drop-in replacement for FastQC | Same output format as FastQC |

For this workshop, we will focus primarily on **FastQC** and **MultiQC**, which together cover the vast majority of short-read QC needs. Long-read QC tools are covered in Section 3.7.

---

## 3.4 FastQC — In Depth

### What is FastQC?

FastQC is a Java-based tool developed by the Babraham Bioinformatics group at the Babraham Institute, Cambridge. It is the most widely used QC tool in NGS and is a standard first step in virtually every short-read analysis pipeline. It reads a FASTQ file and produces an HTML report with a series of graphs and summary statistics, each module checking a different aspect of the data.

### Installation

```bash
# Using conda (recommended)
conda install -c bioconda fastqc

# Using apt (Ubuntu/Debian)
sudo apt install fastqc

# Manual download
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
chmod +x FastQC/fastqc
```

### Basic Usage

```bash
# Run FastQC on a single FASTQ file
fastqc sample_R1.fastq.gz

# Run on multiple files simultaneously
fastqc sample_R1.fastq.gz sample_R2.fastq.gz

# Specify output directory
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o ./fastqc_results/

# Run on all FASTQ files in a directory
fastqc *.fastq.gz -o ./fastqc_results/

# Use multiple threads to process files in parallel
fastqc *.fastq.gz -o ./fastqc_results/ -t 8
```

### Output Files

For each input FASTQ file, FastQC produces two output files:

- `sample_R1_fastqc.html` — the human-readable report; open this in a web browser
- `sample_R1_fastqc.zip` — a zip archive containing the raw data behind each plot; used by MultiQC

### FastQC Summary Icons

Each module in the FastQC report is given one of three status icons:

| Icon | Meaning |
|---|---|
| ✅ Green tick | PASS — results are within normal expected range |
| ⚠️ Orange warning | WARN — results are slightly outside expected range; may or may not be a problem depending on context |
| ❌ Red cross | FAIL — results are well outside expected range; requires investigation |

> ⚠️ **Important:** FastQC's pass/warn/fail thresholds are based on a generic "random library" model. Many legitimate NGS datasets will trigger warnings or failures that are completely expected given the biology of the experiment. Always interpret FastQC results in the context of your experiment type — not blindly by the colour of the icon.

---

## 3.5 Interpreting FastQC Reports — Module by Module

### Module 1: Basic Statistics

This module provides a simple summary of the file:

```
Filename:           sample_R1.fastq.gz
File type:          Conventional base calls
Encoding:           Sanger / Illumina 1.9
Total Sequences:    47,382,951
Sequences flagged as poor quality: 0
Sequence length:    150
%GC:                49
```

**What to check:**
- **Total sequences** — does the number of reads match what you expected from the sequencing order?
- **Sequence length** — is it the read length you requested (e.g. 150 bp)?
- **%GC** — does it roughly match the expected GC content of your organism? Human genome is ~41% GC; if you see 60%, something may be wrong.
- **Sequences flagged as poor quality** — ideally zero; large numbers here indicate a poor sequencing run.

---

### Module 2: Per-Base Sequence Quality

This is the most important plot in the FastQC report. It shows the distribution of Phred quality scores across all reads at each base position.

**What a good plot looks like:**
- The yellow boxes (interquartile range) and whiskers should sit comfortably in the green zone (Q ≥ 28)
- The median line (red) should be ≥ Q30 across most of the read length
- A gradual decline towards the 3' end is normal and expected (see Module 2.6)

**Common problems:**
- Sharp drop in quality at the 3' end — normal; resolved by trimming
- Quality poor from the very start of the read — suggests library or sequencing run problems
- Sudden drop at a specific cycle — may indicate a sequencing run issue (e.g. bubble on the flow cell)
- Consistently low quality throughout — indicative of a failed or very poor sequencing run; contact your sequencing facility

---

### Module 3: Per-Sequence Quality Scores

This plot shows the distribution of mean quality scores across all reads — i.e. for each read, what is its average Q score?

**What a good plot looks like:**
- A single, sharp peak at high quality (Q35–40)
- The vast majority of reads should have a mean Q score above Q30

**Common problems:**
- A bimodal distribution (two peaks) suggests a mixture of high and low quality reads, possibly from separate libraries or a problem with part of the flow cell
- A peak at low Q scores (< Q20) indicates widespread poor base calling

---

### Module 4: Per-Base Sequence Content

This plot shows the proportion of each of the four bases (A, T, G, C) at every position across all reads.

**What a good plot looks like:**
- For random whole-genome libraries, the four lines should run roughly parallel and flat throughout the read, reflecting the genome's base composition
- A and T should be approximately equal; G and C should be approximately equal

**Common problems and expected exceptions:**

- **First 10–15 bases show sequence bias** — this is extremely common and is caused by random hexamer priming bias during cDNA synthesis in RNA-seq, or by ligation bias during library preparation. FastQC will flag this as a FAIL, but it is **expected and generally acceptable** in RNA-seq data. It is usually not trimmed.
- **All reads begin with the same sequence** — indicates amplicon sequencing or a targeted protocol; expected in those contexts
- **%G suddenly increases at a specific position** — may indicate adapter contamination (poly-G tails are common in 2-colour chemistry platforms like NovaSeq when no signal is detected)

---

### Module 5: Per-Sequence GC Content

This plot shows the GC content distribution across all reads, overlaid with the theoretical distribution expected for a random library of the same GC content.

**What a good plot looks like:**
- A smooth, roughly normal (bell-shaped) distribution
- Closely matching the theoretical curve (blue line)

**Common problems:**
- **Sharp spike(s) on top of the main distribution** — indicates overrepresented sequences (e.g. adapter dimers, rRNA contamination, highly expressed transcripts in RNA-seq). These reads all have the same GC content, creating a spike.
- **Broader or shifted distribution** — may indicate contamination with another organism's DNA or systematic library preparation bias
- **Bimodal distribution** — suggests contamination with sequences from a second organism

---

### Module 6: Per-Base N Content

This plots the percentage of bases called as N (undetermined) at each position. The base caller assigns N when it cannot confidently determine which base is present.

**What a good plot looks like:**
- Near-zero N content at all positions (the plot should be essentially flat along the bottom)

**Common problems:**
- Elevated N at the beginning of reads — can occur due to poor calibration of the sequencer at the start of a run; often acceptable if it affects only the first 1–2 bases
- Elevated N throughout — indicates a poor quality sequencing run

---

### Module 7: Sequence Length Distribution

This plot shows the distribution of read lengths across all reads in the file.

**What a good plot looks like:**
- A single sharp peak at the expected read length (e.g. all reads are exactly 150 bp for a 2×150 Illumina run)

**Common problems:**
- Variable read lengths — indicates that some trimming has already been applied to the data, or that the data comes from a platform that naturally produces variable-length reads
- Shorter than expected reads — may indicate quality trimming was applied upstream, or degraded input material

---

### Module 8: Sequence Duplication Levels

This plot shows the proportion of reads that appear multiple times in the dataset. FastQC samples the first 100,000 reads and counts how many times each unique sequence appears.

**What a good plot looks like:**
- For whole genome sequencing, the majority of reads should appear only once (low duplication)
- The blue line (total sequences) and red line (deduplicated sequences) should be close together

**Common problems and expected exceptions:**

- **High duplication in WGS** (> 20–30%) — suggests over-amplification during library preparation, insufficient input DNA, or a problem with the library
- **High duplication in RNA-seq** — this is **normal and expected**. Highly expressed genes generate many identical reads from the same transcripts. Do not be alarmed by high duplication in RNA-seq data.
- **High duplication in targeted/amplicon sequencing** — also expected, as the same amplicon regions are repeatedly sequenced at high depth

> 💡 **Tip:** Duplication rate is better assessed using dedicated tools like **Picard MarkDuplicates** after alignment, where duplicate reads are identified by their mapping coordinates rather than raw sequence, giving a more accurate picture.

---

### Module 9: Overrepresented Sequences

This module lists any sequences that appear in more than 0.1% of all reads — far more than would be expected by chance in a random library. FastQC attempts to identify what these sequences are by comparing them against a built-in database of common contaminants.

**Common sources of overrepresented sequences:**
- **Adapter sequences** — the most common cause; indicates short inserts where the sequencer reads into the adapter
- **rRNA sequences** — common in RNA-seq when rRNA depletion was incomplete
- **PhiX control sequences** — Illumina adds PhiX as a sequencing control; small amounts are expected; large amounts suggest a problem with the run
- **Highly expressed transcripts** — in RNA-seq, a few very highly expressed genes (e.g. globin in blood samples) can dominate

---

### Module 10: Adapter Content

This plot shows the cumulative percentage of reads that contain adapter sequences at each position. FastQC checks for a built-in list of common Illumina adapter sequences.

**What a good plot looks like:**
- Near-zero adapter content throughout the entire read length

**Common problems:**
- **Rising adapter content towards the 3' end** — classic sign of short inserts where the sequencer reads into the adapter. The shorter the average insert size, the earlier the adapter appears and the more reads are affected.
- **Adapters detected from the very first base** — indicates very short inserts or adapter dimers; usually accompanied by high duplication levels

**This module directly informs your trimming strategy in Module 4.** If adapter content is detected, trimming is mandatory before alignment.

---

## 3.6 MultiQC — Aggregating Reports Across Samples

### What is MultiQC?

FastQC is excellent for understanding a single sample, but most NGS projects involve tens to hundreds of samples. Manually opening and comparing 96 individual HTML files is impractical. **MultiQC** solves this problem by automatically scanning a directory for output from FastQC (and over 100 other bioinformatics tools) and compiling everything into a single interactive HTML report.

MultiQC was developed by Phil Ewels at SciLifeLab, Stockholm, and is now one of the most widely cited bioinformatics tools in the field.

### Installation

```bash
# Using conda (recommended)
conda install -c bioconda multiqc

# Using pip
pip install multiqc
```

### Basic Usage

```bash
# Run MultiQC in a directory containing FastQC output
multiqc fastqc_results/

# Specify output directory and report name
multiqc fastqc_results/ -o multiqc_output/ -n my_project_QC_report

# Scan multiple directories
multiqc fastqc_results/ trimming_results/ alignment_results/

# Force overwrite of existing report
multiqc fastqc_results/ -f
```

### What MultiQC Produces

MultiQC generates:
- `multiqc_report.html` — the interactive report; open in a browser
- `multiqc_data/` — a directory of tab-separated files containing the raw data behind each plot, useful for downstream parsing or record keeping

### Key Features of the MultiQC Report

**General Statistics Table**  
At the top of the report is a summary table with one row per sample, showing key metrics side-by-side: total reads, % duplicates, % GC, % adapter content, mean quality score, and more. This allows instant identification of outlier samples.

**Interactive Plots**  
Each FastQC module is presented as a single overlaid plot for all samples. You can click on individual samples to highlight or hide them, zoom in on regions of interest, and export plots for publication.

**Colour-coded Status**  
The status panel at the top shows which FastQC modules passed, warned, or failed for each sample in a compact heatmap, allowing rapid identification of problematic samples before drilling into the details.

### Practical Workflow

```bash
# Step 1: Run FastQC on all samples
fastqc raw_data/*.fastq.gz -o fastqc_results/ -t 16

# Step 2: Aggregate with MultiQC
multiqc fastqc_results/ -o multiqc_output/ -n project_QC

# Step 3: Open the report
open multiqc_output/project_QC.html   # macOS
xdg-open multiqc_output/project_QC.html  # Linux
```

---

## 3.7 QC for Long-Read Data

Long-read data from Oxford Nanopore or PacBio requires different QC tools, as FastQC is designed for short, uniform-length reads and does not handle the variable-length, higher-error-rate characteristics of long-read data well.

### NanoStat

Provides a simple text summary of key statistics for Nanopore FASTQ or FAST5 files.

```bash
# Install
conda install -c bioconda nanostat

# Run on a FASTQ file
NanoStat --fastq nanopore_reads.fastq.gz

# Run on a FAST5 directory
NanoStat --fast5 fast5_directory/ --threads 8
```

**Example NanoStat output:**

```
General summary:
Mean read length:              8,432.5
Mean read quality:                12.3
Median read length:            6,891.0
Median read quality:              13.1
Number of reads:           1,423,891
Read length N50:              14,203
Total bases:           12,001,234,567
Number, percentage and megabases of reads above quality cutoffs
>Q5:   1,423,891 (100.0%) 12,001.2 Mb
>Q7:   1,389,102 (97.6%)  11,734.5 Mb
>Q10:  1,201,043 (84.4%)  10,132.2 Mb
>Q12:    987,234 (69.4%)   8,341.1 Mb
>Q15:    312,891 (22.0%)   2,641.3 Mb
```

**Key metrics for Nanopore data:**
- **Mean/median read quality** — for older R9.4.1 chemistry, Q10–Q14 is typical; for newer R10.4.1 with super accuracy basecalling, Q20+ is achievable
- **Read length N50** — the read length such that 50% of all sequenced bases are in reads of this length or longer; a key measure of library fragment size
- **Total bases** — total yield in gigabases; check against your coverage requirements

### NanoPlot

Produces a rich HTML report with visualisations of read length and quality distributions.

```bash
# Install
conda install -c bioconda nanoplot

# Run on FASTQ
NanoPlot --fastq nanopore_reads.fastq.gz -o nanoplot_output/ --threads 8

# Summary plot only (faster)
NanoPlot --fastq nanopore_reads.fastq.gz -o nanoplot_output/ --plots dot
```

**Key plots in the NanoPlot report:**
- **Read length histogram** — shows the distribution of fragment sizes; should reflect your library preparation (e.g. a size-selected library will have a narrower peak)
- **Read quality histogram** — distribution of per-read mean Q scores
- **Read length vs quality scatter plot** — longer reads often have slightly lower quality; this plot reveals if very long reads are disproportionately poor quality
- **Yield over time** — shows sequencing output over the run; a plateau suggests pore exhaustion or run quality decline

### PacBio QC

For PacBio HiFi (CCS) data, the CCS reads delivered by the sequencing facility have already been through internal quality filtering. Key QC metrics to check:

- **CCS read accuracy** — should be ≥ Q20 (99%) for HiFi reads; typically Q30+ in practice
- **Subread pass number** — the number of times each molecule was sequenced; higher passes = higher accuracy; HiFi typically requires ≥ 3 passes
- **ZMW yield** — the proportion of zero-mode waveguides that produced usable CCS reads

---

## 3.8 Hands-on Exercises

### Exercise 3.1 — Running FastQC on a Single Sample

Download a publicly available test dataset from the NCBI SRA and run FastQC.

```bash
# Create working directory
mkdir -p ~/ngs_workshop/module3 && cd ~/ngs_workshop/module3

# Download a small test FASTQ file (human RNA-seq, SRR subset)
# Option 1: Use SRA tools
conda install -c bioconda sra-tools
fastq-dump --split-files --gzip -X 500000 SRR6821753 -O raw_data/

# Option 2: Download pre-prepared test files (if provided by instructor)
# wget http://[instructor-provided-url]/test_R1.fastq.gz -O raw_data/test_R1.fastq.gz
# wget http://[instructor-provided-url]/test_R2.fastq.gz -O raw_data/test_R2.fastq.gz

# Make output directory
mkdir -p fastqc_results

# Run FastQC
fastqc raw_data/*.fastq.gz -o fastqc_results/ -t 4

# List output files
ls -lh fastqc_results/
```

**Expected output:**
```
test_R1_fastqc.html
test_R1_fastqc.zip
test_R2_fastqc.html
test_R2_fastqc.zip
```

Open `test_R1_fastqc.html` in your browser and work through each module.

---

### Exercise 3.2 — Interpreting a FastQC Report

Work through the following checklist for the FastQC report you just generated:

```
□ Basic Statistics
  □ How many reads are in the file?
  □ What is the read length?
  □ What is the %GC content? Is this expected for human data (~41%)?

□ Per-Base Sequence Quality
  □ Does quality remain above Q30 across most of the read?
  □ At what position does quality begin to drop significantly?
  □ Would you trim this data, and if so, from which position?

□ Per-Sequence GC Content
  □ Does the distribution match the theoretical curve?
  □ Are there any spikes suggesting contamination?

□ Sequence Duplication Levels
  □ What is the estimated duplicate rate?
  □ Is this level of duplication expected for this experiment type?

□ Adapter Content
  □ Is adapter content detected?
  □ From which position does it appear?
  □ Which adapter sequence is present?

□ Overrepresented Sequences
  □ Are any sequences flagged?
  □ What does FastQC identify them as?
```

---

### Exercise 3.3 — Running MultiQC on Multiple Samples

```bash
# If you have multiple FastQC results, run MultiQC to aggregate them
cd ~/ngs_workshop/module3

# Run MultiQC on the FastQC output directory
multiqc fastqc_results/ -o multiqc_output/ -n workshop_QC_report

# Open the report
xdg-open multiqc_output/workshop_QC_report.html
```

**Questions to answer from the MultiQC report:**
1. Which sample has the highest duplicate rate?
2. Which sample has the lowest mean quality score?
3. Are there any samples that stand out as outliers in the general statistics table?
4. Do all samples show adapter contamination, or only some?

---

### Exercise 3.4 — Checking for Contamination with FastQ Screen

```bash
# Install FastQ Screen
conda install -c bioconda fastq-screen

# Download reference databases (this takes time - instructor may provide pre-built databases)
fastq_screen --get_genomes

# Run FastQ Screen
fastq_screen --aligner bowtie2 raw_data/test_R1.fastq.gz --outdir fastqscreen_results/
```

FastQ Screen maps a subset of your reads against a panel of genomes (human, mouse, Ecoli, PhiX, etc.) and reports the percentage that map to each. Ideally, human reads should map almost exclusively to the human genome. High mapping to E. coli may indicate contamination from library preparation reagents.

---

## 3.9 QC Decision Framework — Pass, Warn, or Fail?

After running FastQC and MultiQC, you need to decide what to do next. Use this decision framework:

```
Per-Base Sequence Quality
├── Q30 across most of read → Proceed
├── Quality drops at 3' end only → Trim in Module 4, then proceed
└── Quality poor throughout → Contact sequencing facility; consider resequencing

Adapter Content
├── No adapters detected → Proceed
└── Adapters detected → Trim adapters in Module 4, then proceed

GC Content
├── Matches expected → Proceed
├── Slight deviation → Note and proceed with caution
└── Large spike or completely wrong distribution → Investigate contamination

Sequence Duplication
├── WGS: < 20% → Good
├── WGS: 20–50% → Acceptable; mark duplicates before variant calling
├── WGS: > 50% → Investigate; possible low-input library problem
├── RNA-seq: Any level → Expected; do not remove duplicates
└── Targeted: Any level → Expected; assess on a per-amplicon basis

Per-Base N Content
├── Near zero → Proceed
└── Elevated → Trim or mask N-containing reads

Overrepresented Sequences
├── None detected → Proceed
├── Adapter dimers → Trim
├── rRNA (RNA-seq) → Note poor rRNA depletion; proceed if < 10% of reads
└── Unknown sequence → BLAST to identify; investigate before proceeding
```

---

## 3.10 Common QC Problems and Their Causes

| Problem | Likely Cause | Action |
|---|---|---|
| Low quality throughout all reads | Failed sequencing run; instrument issue | Contact sequencing facility; request rerun |
| Quality poor from cycle 1 | Low cluster density or over-clustering; phasing errors | Contact sequencing facility |
| High adapter content | Short library inserts; size selection failed | Trim adapters (Module 4) |
| Very high duplication (WGS > 60%) | Insufficient input DNA; too many PCR cycles | Flag; consider remaking library with more input |
| GC spike in GC plot | rRNA contamination (RNA-seq); overrepresented sequences | Check overrepresented sequences module; consider rRNA depletion improvement |
| Bimodal quality distribution | Mixed library of two quality populations; flow cell issue | Investigate; may need to split analysis |
| Poly-G tails at read ends | 2-colour chemistry (NovaSeq, NextSeq) with short inserts — no signal reads as G | Trim poly-G tails with fastp or Cutadapt |
| Per-base content failure at read start | Random hexamer priming bias (RNA-seq) or transposase insertion bias (ATAC-seq) | Expected; do not trim unless very severe |
| High %N at specific position | Flow cell bubble or imaging issue at that cycle | Trim to remove affected position |
| Unexpected organism in FastQ Screen | Sample contamination or reagent contamination | Investigate source; filter contaminating reads if necessary |

---

> 💡 **Key takeaway for this module:** FastQC and MultiQC give you a detailed picture of your data's quality, but interpreting that picture requires context. A FAIL in FastQC does not always mean bad data — it means the data differs from a generic expectation. Always ask: *"Given my experiment type, is this expected?"* A high duplication rate in RNA-seq is normal. Adapter content always needs trimming. Per-base content bias at the start of RNA-seq reads is expected and harmless. Develop this critical thinking before you ever run an aligner.

---

**Previous:** [Module 2 — NGS Workflow and Data Generation](./module2_ngs_workflow.md)  
**Next:** [Module 4 — Data Preprocessing and Cleaning](./module4_preprocessing.md)
