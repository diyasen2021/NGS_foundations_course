# Module 4: Data Preprocessing and Cleaning

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 1 | Module 4 of 4 | 🧪 Hands-on

---

## Table of Contents

- [4.1 Why Preprocessing is Necessary](#41-why-preprocessing-is-necessary)
- [4.2 What Does Preprocessing Actually Do?](#42-what-does-preprocessing-actually-do)
- [4.3 Tool Overview](#43-tool-overview)
- [4.4 Trimmomatic — In Depth](#44-trimmomatic--in-depth)
- [4.5 fastp — In Depth](#45-fastp--in-depth)
- [4.6 Cutadapt — In Depth](#46-cutadapt--in-depth)
- [4.7 When NOT to Trim](#47-when-not-to-trim)
- [4.8 Preprocessing for Special Experiment Types](#48-preprocessing-for-special-experiment-types)
- [4.9 Hands-on Exercises](#49-hands-on-exercises)
- [4.10 Evaluating Preprocessing Results](#410-evaluating-preprocessing-results)
- [4.11 Building a Reproducible Preprocessing Pipeline](#411-building-a-reproducible-preprocessing-pipeline)

---

## 4.1 Why Preprocessing is Necessary

In Module 3, you learned how to assess the quality of raw sequencing data. In this module, you act on those findings. Raw FASTQ files are almost never analysis-ready straight out of the sequencer — they contain a variety of contaminants and low-quality data that, if left unaddressed, will degrade every downstream step.

The consequences of skipping preprocessing are concrete and serious:

- **Adapter sequences** in reads will fail to align to the reference genome, reducing mapping rates and potentially introducing alignment artefacts near the true insert boundaries
- **Low-quality bases at read ends** reduce alignment accuracy and can introduce false positive variant calls — a low-quality base miscalled as a SNP looks identical to a real SNP before filtering
- **Poly-G tails** (common on NovaSeq and NextSeq instruments with 2-colour chemistry) will similarly fail to align and inflate soft-clipping rates in BAM files
- **Very short reads** produced after trimming that are below ~30–36 bp often align to multiple locations in the genome and introduce noise into downstream analyses

Preprocessing is a targeted cleaning step. The goal is to remove technical artefacts introduced during library preparation and sequencing, while retaining as much genuine biological signal as possible. It is not about discarding data aggressively — it is about making the data that remains trustworthy.

---

## 4.2 What Does Preprocessing Actually Do?

Modern preprocessing tools can perform several distinct operations, which can be applied individually or in combination:

### Adapter Trimming

Removes adapter sequences that appear at the ends of reads. This happens when the sequenced DNA insert is shorter than the read length, causing the sequencer to read into the adapter ligated to the other end of the fragment.

```
Read (150 bp):  [----insert (80 bp)----][----adapter (70 bp)----]
After trimming: [----insert (80 bp)----]
```

Most tools can detect adapters automatically by looking for sequences that appear frequently at read ends, or you can supply the adapter sequence explicitly.

### Quality Trimming

Removes bases from the ends of reads where the Phred quality score drops below a defined threshold. The two main approaches are:

- **Hard cutoff** — trim all bases below Q20 from the 3' end
- **Sliding window** — scan along the read in a window (e.g. 4 bases); trim from the point where the average quality in the window drops below a threshold (e.g. Q20). This is more robust than hard cutoffs because a single good base surrounded by poor bases will not rescue a bad region.

### Minimum Length Filtering

After trimming, some reads become very short. Reads shorter than a defined minimum (typically 30–36 bp) are discarded entirely, as they are too short to align reliably and uniquely to the genome.

### Quality Filtering

Discard entire reads (or read pairs) where the overall quality is too low — for example, reads where more than 40% of bases are below Q20, or reads with a mean quality below Q15.

### Poly-X Tail Trimming

Removes stretches of the same base at read ends — most commonly poly-G tails on NovaSeq/NextSeq instruments, where a lack of signal (no incorporation) is interpreted as G by the 2-colour chemistry. Also removes poly-A tails that appear in RNA-seq reads when the sequencer reads through the poly-A tail of a transcript.

### N Base Filtering

Discards reads containing more than a threshold number of uncalled bases (N), as these positions contribute no reliable information and can cause alignment problems.

### Paired-End Coordination

For paired-end data, if one read in a pair is discarded (e.g. it becomes too short after trimming), the other read must also be handled — either discarded (strict pairing) or written to an "unpaired" output file. Maintaining proper pairing is critical because most aligners and downstream tools expect perfectly matched pairs.

---

## 4.3 Tool Overview

| Tool | Language | Speed | Best For | Auto-detect Adapters |
|---|---|---|---|---|
| **Trimmomatic** | Java | Moderate | Illumina PE/SE; highly configurable | No — must specify |
| **fastp** | C++ | Very fast | General purpose; great default settings | Yes |
| **Cutadapt** | Python | Fast | Flexible adapter trimming; miRNA-seq; custom workflows | Partial |
| **TrimGalore** | Perl wrapper | Moderate | RNA-seq; wraps Cutadapt with sensible defaults | Yes |
| **BBDuk** | Java (BBTools) | Very fast | Contamination removal + trimming combined | Yes |

For this workshop we focus on **Trimmomatic**, **fastp**, and **Cutadapt** — the three most widely encountered tools in Indian genomics labs and published pipelines. fastp is increasingly the recommended choice for new projects due to its speed, automatic adapter detection, and excellent HTML reports.

---

## 4.4 Trimmomatic — In Depth

### Overview

Trimmomatic was developed at the Usadel Lab (RWTH Aachen) and published in 2014. Despite being written in Java and requiring you to specify adapters explicitly, it remains one of the most widely cited preprocessing tools and is used in many established pipelines including nf-core/rnaseq (as an option) and many GATK-based workflows.

Trimmomatic processes reads in a sequential pipeline of steps, applying each operation in the order you specify on the command line. The order matters — adapter trimming should always come before quality trimming.

### Installation

```bash
# Using conda (recommended)
conda install -c bioconda trimmomatic

# Manual download
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
```

### Adapter Files

Trimmomatic requires you to provide the adapter sequences in a FASTA file. Standard adapter files are bundled with the Trimmomatic installation:

```bash
# Find the adapter files after conda install
ls $(dirname $(which trimmomatic))/../share/trimmomatic/adapters/

# Common adapter files:
# TruSeq2-PE.fa       — older Illumina paired-end (HiSeq 2000, GAII)
# TruSeq3-PE.fa       — newer Illumina paired-end (HiSeq 2500, MiSeq, NovaSeq)
# TruSeq3-PE-2.fa     — TruSeq3 with additional sequences
# TruSeq3-SE.fa       — single-end TruSeq3
# NexteraPE-PE.fa     — Nextera XT paired-end
```

If you are unsure which adapters were used, check with your sequencing facility, or use a tool like fastp (which auto-detects) to confirm first.

### Paired-End Usage (Most Common)

```bash
trimmomatic PE \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_paired.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_paired.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36
```

### Single-End Usage

```bash
trimmomatic SE \
    sample.fastq.gz \
    sample_trimmed.fastq.gz \
    ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36
```

### Understanding the Parameters

**`ILLUMINACLIP:adapter_file:seed_mismatches:palindrome_clip_threshold:simple_clip_threshold`**

| Parameter | Value | Meaning |
|---|---|---|
| `adapter_file` | TruSeq3-PE-2.fa | Path to FASTA file containing adapter sequences |
| `seed_mismatches` | 2 | Maximum mismatches allowed in the 16 bp adapter seed match |
| `palindrome_clip_threshold` | 30 | Score threshold for palindrome adapter detection (PE mode); 30 is recommended |
| `simple_clip_threshold` | 10 | Score threshold for simple adapter detection; 10 is recommended |

**`LEADING:3`** — Remove bases from the start (5' end) of a read if their quality is below 3. Catches any very poor bases at read start.

**`TRAILING:3`** — Remove bases from the end (3' end) of a read if their quality is below 3.

**`SLIDINGWINDOW:4:20`** — Scan with a 4-base sliding window; trim when the average quality within the window drops below 20. Applied after LEADING/TRAILING.

**`MINLEN:36`** — Discard any read shorter than 36 bases after all trimming steps.

**`HEADCROP:N`** — Remove the first N bases from every read regardless of quality. Useful for removing fixed-length primer sequences or correcting per-base content bias at the start of reads.

**`CROP:N`** — Trim reads to a maximum length of N bases. Useful when you want all reads to be the same length.

### Trimmomatic Output Summary

Trimmomatic prints a summary to the terminal after completion:

```
TrimmomaticPE: Started with arguments:
Input Read Pairs: 47382951
Both Surviving: 45891203 (96.86%)
Forward Only Surviving: 891042 (1.88%)
Reverse Only Surviving: 312847 (0.66%)
Dropped: 287859 (0.61%)
TrimmomaticPE: Completed successfully
```

**What to look for:**
- **Both Surviving** should ideally be > 90%; lower values suggest aggressive trimming settings or poor data quality
- **Dropped** rates above 5–10% warrant investigation — check whether your trimming parameters are too aggressive

### Running Trimmomatic with Multiple Threads

```bash
trimmomatic PE -threads 8 \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_paired.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_paired.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
```

---

## 4.5 fastp — In Depth

### Overview

fastp was developed by Shifu Chen at HaploX Genomics and published in 2018. It is written in C++ and is dramatically faster than Trimmomatic — typically 3–5× faster on the same data. Its key advantages over Trimmomatic are:

- **Automatic adapter detection** — no need to specify adapter files; fastp detects adapters from the data itself by analysing the overlap between paired reads
- **Built-in HTML and JSON reports** — generates quality reports for both pre- and post-trimming data in a single run, removing the need to run FastQC separately (though running FastQC after trimming is still good practice)
- **Poly-G and poly-X trimming** — handles NovaSeq/NextSeq poly-G tails automatically
- **Simpler command line** — sensible defaults mean you often need very few arguments

fastp is now the recommended choice for most new short-read preprocessing workflows.

### Installation

```bash
# Using conda (recommended)
conda install -c bioconda fastp

# Using apt (Ubuntu)
sudo apt install fastp
```

### Basic Paired-End Usage

```bash
# Minimal command — fastp auto-detects adapters
fastp \
    -i sample_R1.fastq.gz \
    -I sample_R2.fastq.gz \
    -o sample_R1_trimmed.fastq.gz \
    -O sample_R2_trimmed.fastq.gz
```

### Full Recommended Usage

```bash
fastp \
    -i sample_R1.fastq.gz \
    -I sample_R2.fastq.gz \
    -o sample_R1_trimmed.fastq.gz \
    -O sample_R2_trimmed.fastq.gz \
    --json fastp_results/sample_fastp.json \
    --html fastp_results/sample_fastp.html \
    --thread 8 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 36 \
    --detect_adapter_for_pe \
    --correction \
    --trim_poly_g \
    --trim_poly_x \
    --report_title "Sample QC Report"
```

### Key fastp Parameters Explained

| Parameter | Default | Meaning |
|---|---|---|
| `--qualified_quality_phred 20` | 15 | Minimum quality score to consider a base "qualified" |
| `--unqualified_percent_limit 40` | 40 | Discard read if > 40% of bases are below the quality threshold |
| `--length_required 36` | 15 | Discard reads shorter than 36 bp after trimming |
| `--detect_adapter_for_pe` | off | Enable automatic adapter detection for paired-end data (recommended) |
| `--adapter_sequence` | auto | Manually specify adapter sequence if known |
| `--correction` | off | Correct mismatched bases in overlapping paired-end reads |
| `--trim_poly_g` | auto | Trim poly-G tails (auto-enabled for NovaSeq/NextSeq data) |
| `--trim_poly_x` | off | Trim any poly-X tail |
| `--thread` | 3 | Number of CPU threads |
| `--json` | fastp.json | Path for JSON report output |
| `--html` | fastp.html | Path for HTML report output |

### fastp HTML Report

The fastp HTML report contains before-and-after comparisons for:
- Per-base quality distribution
- Read quality distribution  
- GC content distribution
- Read length distribution after trimming
- Adapter content detected and removed
- Filtering statistics (reads passed, failed, with adapter, with low quality, too short)

This makes fastp particularly powerful for a workshop setting — one command produces both the trimmed output and the QC report.

### fastp Summary Statistics

```
Read1 before filtering:
total reads: 47382951
total bases: 7107442650
Q20 bases: 6931047291 (97.52%)
Q30 bases: 6612241567 (93.03%)

Read1 after filtering:
total reads: 46201834
total bases: 6817948123
Q20 bases: 6798231012 (99.71%)
Q30 bases: 6593214891 (96.71%)

Filtering result:
reads passed filter: 46201834 (97.51%)
reads failed due to low quality: 892341 (1.88%)
reads failed due to too many N: 12093 (0.03%)
reads failed due to too short: 276683 (0.58%)
reads with adapter trimmed: 8234921 (17.38%)
bases trimmed due to adapters: 142934821
```

---

## 4.6 Cutadapt — In Depth

### Overview

Cutadapt was developed by Marcel Martin and is one of the oldest and most flexible adapter trimming tools available. While fastp and Trimmomatic are better suited to general WGS and RNA-seq workflows, Cutadapt excels in situations requiring precise, custom adapter trimming:

- **Small RNA sequencing (miRNA-seq)** — requires very specific 3' adapter trimming and handling of reads shorter than the adapter itself
- **Custom library protocols** with non-standard adapter sequences
- **CLIP-seq, RIP-seq** and other specialised protocols
- **Situations where you need fine-grained control** over exactly what is trimmed and how

### Installation

```bash
# Using conda
conda install -c bioconda cutadapt

# Using pip
pip install cutadapt
```

### Basic Usage

```bash
# Single-end, trim 3' adapter
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -o sample_trimmed.fastq.gz \
    sample.fastq.gz

# Paired-end, trim 3' adapters from both reads
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o sample_R1_trimmed.fastq.gz \
    -p sample_R2_trimmed.fastq.gz \
    sample_R1.fastq.gz sample_R2.fastq.gz

# With quality trimming and minimum length
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20 \
    --minimum-length 36 \
    -o sample_R1_trimmed.fastq.gz \
    -p sample_R2_trimmed.fastq.gz \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    --cores 8
```

### Common Illumina Adapter Sequences

```
# Illumina TruSeq Universal Adapter (Read 1 3' adapter)
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# Illumina TruSeq Read 2 3' adapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# Nextera/Tn5 Read 1 adapter
CTGTCTCTTATACACATCT

# Nextera/Tn5 Read 2 adapter
CTGTCTCTTATACACATCT

# Small RNA 3' adapter (TruSeq Small RNA)
TGGAATTCTCGGGTGCCAAGG
```

### Cutadapt for miRNA-seq

Small RNA sequencing is a special case where the reads are deliberately short (miRNAs are 18–25 nt), meaning reads almost always contain adapter sequence. Cutadapt handles this well:

```bash
cutadapt \
    -a TGGAATTCTCGGGTGCCAAGG \
    --discard-untrimmed \
    --minimum-length 16 \
    --maximum-length 35 \
    -q 20 \
    -o mirna_trimmed.fastq.gz \
    mirna_raw.fastq.gz \
    --cores 8
```

`--discard-untrimmed` removes reads where no adapter was found — in miRNA-seq, a read without an adapter is likely noise or contamination rather than a genuine small RNA.

---

## 4.7 When NOT to Trim

Trimming is not always the right choice. There are situations where trimming can actually harm your analysis:

### Do Not Aggressively Quality-Trim Before Alignment with Soft-Clipping Aligners

Modern aligners like **STAR** and **HISAT2** use soft-clipping, meaning they can handle low-quality ends by ignoring them during alignment. Over-trimming RNA-seq reads can reduce the number of reads that span exon-exon junctions, making splice junction detection less sensitive. For RNA-seq, light or no quality trimming is often recommended when using these aligners — adapter trimming is still necessary.

### Do Not Remove Duplicates from RNA-seq Data

This is not a trimming step, but it bears repeating here: do **not** run deduplication on RNA-seq data. Duplicate reads in RNA-seq are genuine — multiple molecules of the same transcript exist in the cell and generate identical reads. Removing them would reduce expression estimates for highly expressed genes and distort differential expression analysis.

### Do Not Trim Per-Base Content Bias at Read Starts

As discussed in Module 3, the per-base sequence content bias at the first 10–15 bases of RNA-seq reads (caused by random hexamer priming) is a technical artefact, but trimming it off removes real sequence data. Modern aligners handle this bias without needing to remove the bases.

### Minimal Trimming for Long Reads

For Oxford Nanopore and PacBio long reads, aggressive quality trimming is counterproductive. Long-read aligners (minimap2) are designed to handle the higher error rates of long reads. Trimming too aggressively shortens reads unnecessarily, eliminating the very length advantage that makes long-read sequencing valuable. Adapter removal is still recommended for Nanopore data (using **Porechop**), but quality trimming is typically not applied.

---

## 4.8 Preprocessing for Special Experiment Types

### RNA-seq

**Recommended approach:** Adapter trimming + minimum length filter. Light quality trimming (SLIDINGWINDOW:4:15 or none).

```bash
fastp \
    -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o sample_R1_trimmed.fastq.gz -O sample_R2_trimmed.fastq.gz \
    --detect_adapter_for_pe \
    --length_required 30 \
    --thread 8 \
    --html fastp_rnaseq.html --json fastp_rnaseq.json
```

### ATAC-seq

ATAC-seq uses the Tn5 transposase, which inserts sequencing adapters directly into open chromatin. Reads often have Nextera adapters and can be very short (nucleosome-free regions are ~150 bp). Trim Nextera adapters aggressively.

```bash
fastp \
    -i atac_R1.fastq.gz -I atac_R2.fastq.gz \
    -o atac_R1_trimmed.fastq.gz -O atac_R2_trimmed.fastq.gz \
    --adapter_sequence CTGTCTCTTATACACATCT \
    --adapter_sequence_r2 CTGTCTCTTATACACATCT \
    --length_required 25 \
    --thread 8
```

### ChIP-seq

Similar to WGS — trim adapters and low quality bases. Short reads are acceptable as most peaks are in mappable regions.

```bash
trimmomatic SE \
    chipseq_IP.fastq.gz chipseq_IP_trimmed.fastq.gz \
    ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
    SLIDINGWINDOW:4:20 MINLEN:25
```

### Bisulfite Sequencing (WGBS / RRBS)

Bisulfite-converted DNA has reduced sequence complexity (most cytosines are converted to thymines), making adapter contamination harder to detect. **TrimGalore** is the standard tool for bisulfite data as it handles the reduced complexity correctly.

```bash
# Install TrimGalore
conda install -c bioconda trim-galore

# Run TrimGalore for RRBS
trim_galore --paired --rrbs \
    bisulfite_R1.fastq.gz bisulfite_R2.fastq.gz \
    -o trimgalore_output/
```

### Nanopore — Adapter Removal with Porechop

```bash
# Install Porechop
conda install -c bioconda porechop

# Remove adapters from Nanopore reads
porechop \
    -i nanopore_reads.fastq.gz \
    -o nanopore_trimmed.fastq.gz \
    --threads 8
```

---

## 4.9 Hands-on Exercises

### Setup

```bash
# Create directory structure for this module
mkdir -p ~/ngs_workshop/module4/{raw_data,trimmed,fastqc_before,fastqc_after,reports}
cd ~/ngs_workshop/module4

# If you completed Module 3, copy your raw data here
cp ~/ngs_workshop/module3/raw_data/*.fastq.gz raw_data/

# Install tools if not already installed
conda install -c bioconda trimmomatic fastp cutadapt fastqc multiqc -y
```

---

### Exercise 4.1 — Run FastQC Before Trimming

Always run FastQC on the raw data first so you have a baseline to compare against after trimming.

```bash
fastqc raw_data/*.fastq.gz -o fastqc_before/ -t 4
```

Note the following from the pre-trimming report before proceeding:
- Adapter content: present / absent, and from which position?
- Per-base quality: where does it drop below Q28?
- Duplication rate
- Any overrepresented sequences?

---

### Exercise 4.2 — Trim with Trimmomatic

```bash
# Find the adapter file path
ADAPTER_PATH=$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-PE-2.fa

# Run Trimmomatic PE
trimmomatic PE -threads 4 \
    raw_data/SRR6821753_1.fastq.gz \
    raw_data/SRR6821753_2.fastq.gz \
    trimmed/SRR6821753_R1_paired.fastq.gz \
    trimmed/SRR6821753_R1_unpaired.fastq.gz \
    trimmed/SRR6821753_R2_paired.fastq.gz \
    trimmed/SRR6821753_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTER_PATH}:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36

# Check output file sizes
ls -lh trimmed/
```

**Questions:**
1. What percentage of read pairs survived trimming?
2. How many reads were dropped entirely?
3. How large are the paired output files compared to the input files?

---

### Exercise 4.3 — Trim with fastp

```bash
# Run fastp on the same raw data
fastp \
    -i raw_data/SRR6821753_1.fastq.gz \
    -I raw_data/SRR6821753_2.fastq.gz \
    -o trimmed/SRR6821753_R1_fastp.fastq.gz \
    -O trimmed/SRR6821753_R2_fastp.fastq.gz \
    --json reports/SRR6821753_fastp.json \
    --html reports/SRR6821753_fastp.html \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 36 \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --thread 4

# Open the fastp HTML report
xdg-open reports/SRR6821753_fastp.html
```

**Questions:**
1. What adapters did fastp detect automatically?
2. What percentage of reads had adapters trimmed?
3. How does the Q30 rate before and after trimming compare?

---

### Exercise 4.4 — Run FastQC After Trimming

```bash
# Run FastQC on trimmed output from both tools
fastqc \
    trimmed/SRR6821753_R1_paired.fastq.gz \
    trimmed/SRR6821753_R2_paired.fastq.gz \
    trimmed/SRR6821753_R1_fastp.fastq.gz \
    trimmed/SRR6821753_R2_fastp.fastq.gz \
    -o fastqc_after/ -t 4
```

---

### Exercise 4.5 — Compare Before and After with MultiQC

```bash
# Run MultiQC across both pre- and post-trimming FastQC results
multiqc fastqc_before/ fastqc_after/ reports/ \
    -o reports/multiqc_comparison \
    -n preprocessing_comparison

xdg-open reports/multiqc_comparison/preprocessing_comparison.html
```

Use the MultiQC report to answer the following:

```
Before vs After Trimming Comparison Checklist:

□ Adapter Content
  □ Was adapter content present before trimming?
  □ Is it gone (or reduced to near zero) after trimming?

□ Per-Base Sequence Quality
  □ Has the quality profile improved at read ends?
  □ Are there any new quality issues introduced by trimming?

□ Sequence Length Distribution
  □ Before trimming: all reads the same length?
  □ After trimming: is there now a distribution of lengths?
  □ What is the new minimum read length?

□ Trimmomatic vs fastp comparison
  □ Which tool retained more reads?
  □ Which tool produced a better quality profile?
  □ Do the two tools agree on adapter sequences?
```

---

### Exercise 4.6 — Write a Preprocessing Shell Script

A good bioinformatician does not run commands manually for each sample. Write a shell script that processes all samples in a directory automatically.

```bash
# Create the script
cat > ~/ngs_workshop/module4/run_preprocessing.sh << 'EOF'
#!/bin/bash
# NGS Preprocessing Pipeline
# Usage: bash run_preprocessing.sh <sample_name> <R1.fastq.gz> <R2.fastq.gz>

set -euo pipefail   # Exit on error, undefined variable, or pipe failure

SAMPLE=$1
R1=$2
R2=$3

OUTDIR="trimmed"
REPORTDIR="reports"
THREADS=4

mkdir -p ${OUTDIR} ${REPORTDIR}

echo "======================================"
echo "Processing sample: ${SAMPLE}"
echo "======================================"

fastp \
    -i ${R1} \
    -I ${R2} \
    -o ${OUTDIR}/${SAMPLE}_R1_trimmed.fastq.gz \
    -O ${OUTDIR}/${SAMPLE}_R2_trimmed.fastq.gz \
    --json ${REPORTDIR}/${SAMPLE}_fastp.json \
    --html ${REPORTDIR}/${SAMPLE}_fastp.html \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 36 \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --thread ${THREADS} \
    --report_title "${SAMPLE} Preprocessing Report"

echo "Done: ${SAMPLE}"
echo "Output: ${OUTDIR}/${SAMPLE}_R1_trimmed.fastq.gz"
echo "Report: ${REPORTDIR}/${SAMPLE}_fastp.html"
EOF

chmod +x ~/ngs_workshop/module4/run_preprocessing.sh

# Run it
bash ~/ngs_workshop/module4/run_preprocessing.sh \
    SRR6821753 \
    raw_data/SRR6821753_1.fastq.gz \
    raw_data/SRR6821753_2.fastq.gz
```

---

## 4.10 Evaluating Preprocessing Results

After trimming, always verify the output before proceeding to alignment. Use these checks:

### Check 1 — File Integrity

```bash
# Verify gzip files are not corrupted
gzip -t trimmed/SRR6821753_R1_trimmed.fastq.gz && echo "R1: OK"
gzip -t trimmed/SRR6821753_R2_trimmed.fastq.gz && echo "R2: OK"
```

### Check 2 — Read Count Consistency

Paired-end files must have exactly the same number of reads after trimming.

```bash
# Count reads in R1
echo "R1 read count:"
zcat trimmed/SRR6821753_R1_trimmed.fastq.gz | wc -l | awk '{print $1/4}'

# Count reads in R2
echo "R2 read count:"
zcat trimmed/SRR6821753_R2_trimmed.fastq.gz | wc -l | awk '{print $1/4}'

# They must match exactly — if they don't, the pairing is broken
```

### Check 3 — FastQC on Trimmed Output

```bash
fastqc trimmed/*_trimmed.fastq.gz -o fastqc_after/ -t 4
```

**What to confirm in the post-trimming FastQC report:**
- Adapter Content module should show near-zero adapter content
- Per-Base Sequence Quality should show improvement at read ends
- Sequence Length Distribution should now show a range of lengths rather than a single value
- All other modules should be the same or better than before trimming

### Check 4 — Data Retention Rate

```bash
# A simple summary using fastp JSON output
python3 -c "
import json
with open('reports/SRR6821753_fastp.json') as f:
    data = json.load(f)
summary = data['summary']
print(f'Input reads: {summary[\"before_filtering\"][\"total_reads\"]:,}')
print(f'Output reads: {summary[\"after_filtering\"][\"total_reads\"]:,}')
retention = summary['after_filtering']['total_reads'] / summary['before_filtering']['total_reads'] * 100
print(f'Retention rate: {retention:.1f}%')
print(f'Q30 before: {summary[\"before_filtering\"][\"q30_rate\"]*100:.1f}%')
print(f'Q30 after: {summary[\"after_filtering\"][\"q30_rate\"]*100:.1f}%')
"
```

---

## 4.11 Building a Reproducible Preprocessing Pipeline

One of the core principles of good bioinformatics is **reproducibility** — another researcher (or your future self) should be able to re-run your analysis and get identical results. This requires careful attention to how you document and structure your preprocessing steps.

### Record Your Tool Versions

```bash
# Always record the exact version of every tool used
fastqc --version
trimmomatic -version
fastp --version
cutadapt --version

# Save to a file for your lab notebook / methods section
echo "Tool versions for preprocessing:" > reports/tool_versions.txt
fastqc --version >> reports/tool_versions.txt
trimmomatic -version >> reports/tool_versions.txt
fastp --version >> reports/tool_versions.txt
```

### Use Conda Environments

```bash
# Create a dedicated environment for your project
conda create -n ngs_workshop python=3.10
conda activate ngs_workshop
conda install -c bioconda fastqc trimmomatic fastp cutadapt multiqc

# Export the environment specification
conda env export > environment.yml

# Anyone can recreate your exact environment with:
# conda env create -f environment.yml
```

### Keep Raw Data Untouched

Never overwrite or modify your original raw FASTQ files. Always write trimmed output to a new directory. If something goes wrong with preprocessing, you need to be able to start over.

```bash
# Recommended directory structure
project/
├── raw_data/           ← original FASTQ files; read-only; never modified
├── trimmed/            ← output of preprocessing
├── aligned/            ← output of alignment (Module 5)
├── variants/           ← output of variant calling (Module 6)
├── reports/            ← all QC reports (FastQC, MultiQC, fastp HTML)
├── logs/               ← log files from each tool
├── scripts/            ← all shell scripts used in the analysis
└── environment.yml     ← conda environment specification
```

### Write Everything to a Log File

```bash
# Redirect both stdout and stderr to a log file
fastp \
    -i raw_data/sample_R1.fastq.gz \
    -I raw_data/sample_R2.fastq.gz \
    -o trimmed/sample_R1_trimmed.fastq.gz \
    -O trimmed/sample_R2_trimmed.fastq.gz \
    --detect_adapter_for_pe \
    --thread 8 \
    --html reports/sample_fastp.html \
    --json reports/sample_fastp.json \
    2> logs/sample_fastp.log

echo "Exit code: $?" >> logs/sample_fastp.log
```

---

> 💡 **Key takeaway for this module:** Preprocessing is a bridge between raw data and meaningful analysis. The goal is not to discard as much data as possible — it is to remove only what is genuinely harmful (adapters, very low quality bases) while preserving as much biological signal as possible. fastp is an excellent default choice for most Illumina short-read experiments. Always run FastQC before and after preprocessing to confirm the cleaning has worked as expected, and always verify that paired-end files contain matched read counts before proceeding to alignment.

---

**Previous:** [Module 3 — Quality Control of Raw Sequencing Data](./module3_quality_control.md)  
**Next:** [Module 5 — Sequence Alignment and Mapping](./module5_alignment.md) *(Day 2)*
