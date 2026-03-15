# Module 5: Sequence Alignment and Mapping

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 2 | Module 5 of 4 | 🧪 Hands-on

> **Context:** This module focuses on DNA alignment for **variant calling**. The goal is to produce a clean, sorted, indexed BAM file that is ready for GATK variant calling in Module 6.

---

## Table of Contents

- [5.1 What is Sequence Alignment?](#51-what-is-sequence-alignment)
- [5.2 The Reference Genome](#52-the-reference-genome)
- [5.3 How BWA-MEM2 Works](#53-how-bwa-mem2-works)
- [5.4 The SAM and BAM File Formats](#54-the-sam-and-bam-file-formats)
- [5.5 Post-Alignment Processing](#55-post-alignment-processing)
- [5.6 Marking Duplicate Reads](#56-marking-duplicate-reads)
- [5.7 Base Quality Score Recalibration (BQSR)](#57-base-quality-score-recalibration-bqsr)
- [5.8 Alignment Quality Metrics](#58-alignment-quality-metrics)
- [5.9 The Full Alignment Pipeline — Summary](#59-the-full-alignment-pipeline--summary)
- [5.10 Hands-on Exercises](#510-hands-on-exercises)

---

## 5.1 What is Sequence Alignment?

After preprocessing in Module 4, you have clean FASTQ files containing millions of short reads — each one a 150 bp snippet of DNA from somewhere in the genome. On their own these reads have no genomic context. You do not know which chromosome they came from, which gene they overlap, or whether a base difference from the expected sequence represents a real variant or a sequencing error.

**Sequence alignment** is the process of mapping each read back to its position of origin in a **reference genome**. Once reads are aligned, every position in the genome has a stack of reads covering it — this is your **coverage**. By examining what bases the reads carry at each position, you can identify locations where the sample differs from the reference — these are your **variants**.

The output of alignment is a **BAM file** — a compressed, binary version of the aligned reads — which becomes the primary input for all downstream variant calling, coverage analysis, and visualisation.

```
FASTQ reads (unplaced)          Reference genome
                                Chr1: ATCGGCTAGCTAGCTAGCT...
Read 1: ATCGGCTAGCT     →       |||||||||||
Read 2:    GGCTAGCTAG   →          ||||||||||||
Read 3:       AGCTAGCTA →             |||||||||
                                      ↓
                                BAM file — reads placed on genome
                                Ready for variant calling
```

---

## 5.2 The Reference Genome

### What is a Reference Genome?

The reference genome is the standard DNA sequence of a species against which your sample reads are aligned. For human genomics, the current standard is **GRCh38** (also called hg38), released by the Genome Reference Consortium. It is important to understand what the reference genome is — and is not:

- It is **not** the genome of a single person. It is a mosaic assembled from the DNA of multiple anonymous donors.
- It is **not** perfect. It contains gaps, regions of low confidence, and known errors that are progressively corrected in each release.
- It is **not** the only valid reference. The newer **T2T-CHM13** (Telomere-to-Telomere) assembly released in 2022 is a more complete representation of the human genome, including previously unsequenced centromeric and telomeric regions.

### Why the Reference Choice Matters

Aligning to GRCh37/hg19 versus GRCh38/hg38 will give you different genomic coordinates for the same variants. **Never mix alignments from different reference versions in the same analysis.** Always confirm which reference was used when comparing your results to published datasets or clinical databases.

### Downloading the Reference Genome

```bash
# Create reference directory
mkdir -p ~/ngs_workshop/reference && cd ~/ngs_workshop/reference

# Download GRCh38 reference (full genome — ~900 MB compressed)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# For the workshop — use chromosome 20 only (much smaller, ~65 MB)
# This lets us run a complete pipeline quickly without waiting hours
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
gunzip chr20.fa.gz
```

> 💡 **Workshop note:** For hands-on exercises we use **chromosome 20 only**. It is large enough to demonstrate all pipeline steps realistically but small enough to process in minutes rather than hours. All commands in this module work identically on the full genome.

### Indexing the Reference

Before alignment, the reference must be **indexed** — a process that creates lookup tables allowing BWA-MEM2 to find candidate alignment positions rapidly without scanning the entire genome for every read.

```bash
# Index the reference with BWA-MEM2 (do this once; takes ~5-10 minutes for chr20)
bwa-mem2 index chr20.fa

# Also create a FASTA index (used by SAMtools and GATK)
samtools faidx chr20.fa

# Create a sequence dictionary (required by GATK)
gatk CreateSequenceDictionary -R chr20.fa
```

This creates several index files alongside the reference:
```
chr20.fa
chr20.fa.fai          ← SAMtools FASTA index
chr20.dict            ← GATK sequence dictionary
chr20.fa.0123         ← BWA-MEM2 index files
chr20.fa.amb
chr20.fa.ann
chr20.fa.bwt.2bit.64
chr20.fa.pac
```

---

## 5.3 How BWA-MEM2 Works

### Why BWA-MEM2 for Variant Calling?

For DNA variant calling, **BWA-MEM2** is the industry standard aligner. It is the recommended aligner in the GATK Best Practices pipeline for germline and somatic variant calling. It is an optimised, faster reimplementation of the original BWA-MEM algorithm.

Key properties that make it suitable for variant calling:
- Designed for **short reads (70 bp – 1 Mb)** aligned to a large reference genome
- Handles **paired-end reads** natively, using the insert size distribution to improve mapping
- Uses **local alignment** — only the best-matching portion of a read is aligned, rather than forcing the entire read to match; this handles soft-clipped ends gracefully
- Produces **mapping quality scores (MAPQ)** that reflect confidence in each alignment — critical for filtering in variant calling

BWA-MEM2 is **not** suitable for RNA-seq (it does not handle splicing) or long reads (use minimap2 for those).

### The Seed-and-Extend Algorithm

BWA-MEM2 uses a two-step approach:

**Step 1 — Seeding:** The read is broken into short exact-match seeds (~19 bp) which are looked up in the FM-index of the reference genome. This rapidly identifies candidate alignment regions without checking every position.

**Step 2 — Extension:** The full read is aligned to each candidate region using Smith-Waterman dynamic programming, allowing mismatches and gaps (insertions and deletions). The best alignment is chosen and assigned a mapping quality score.

### Mapping Quality (MAPQ)

Every aligned read receives a MAPQ score — a Phred-scaled probability that the read is mapped to the wrong location:

| MAPQ | Meaning |
|---|---|
| 0 | Read maps equally well to multiple locations — position uncertain |
| 1–19 | Low confidence alignment |
| 20–59 | Moderate to high confidence |
| 60 | Maximum confidence — read maps uniquely and perfectly |

For variant calling, reads with MAPQ < 20 are typically filtered out. Reads mapping to repetitive regions (SINEs, LINEs, centromeres) often have MAPQ = 0 and contribute nothing to variant calls.

---

## 5.4 The SAM and BAM File Formats

### SAM Format

SAM (Sequence Alignment/Map) is a tab-delimited text format for storing aligned reads. Every aligned read occupies one line with 11 mandatory fields:

```
QNAME   FLAG  RNAME  POS    MAPQ  CIGAR   RNEXT  PNEXT  TLEN  SEQ         QUAL
read1   99    chr20  61839  60    150M    =      62012  323   ATCG...     BBFF...
```

| Field | Name | Description |
|---|---|---|
| QNAME | Query name | The read name from the FASTQ file |
| FLAG | Bitwise flag | Encodes read properties (paired, mapped, reverse strand, etc.) |
| RNAME | Reference name | Chromosome the read mapped to (e.g. chr20) |
| POS | Position | 1-based leftmost mapping position on the chromosome |
| MAPQ | Mapping quality | Phred-scaled confidence in the alignment |
| CIGAR | CIGAR string | Describes the alignment (matches, insertions, deletions, clips) |
| SEQ | Sequence | The nucleotide sequence of the read |
| QUAL | Quality | Per-base Phred quality scores (same as FASTQ line 4) |

### Understanding CIGAR Strings

The CIGAR string describes exactly how the read aligns to the reference:

| Code | Meaning | Example |
|---|---|---|
| M | Alignment match (can be match or mismatch) | 150M = 150 bp aligned |
| I | Insertion in read relative to reference | 5I = 5 bp inserted |
| D | Deletion in read relative to reference | 3D = 3 bp deleted |
| S | Soft clip — bases present in read but not aligned | 10S140M = first 10 bp clipped |
| N | Skipped region (used in RNA-seq for introns) | Not used in DNA alignment |

Examples:
```
150M          → perfect 150 bp alignment, no gaps
130M5I15M     → 130 bp match, 5 bp insertion, 15 bp match
140M10S       → 140 bp aligned, last 10 bp soft-clipped (low quality or adapter)
```

### Understanding FLAG Values

The FLAG field is a single integer that encodes multiple properties as binary bits. Common FLAG values:

| FLAG | Meaning |
|---|---|
| 0 | Read mapped, forward strand, single-end |
| 16 | Read mapped, reverse strand |
| 99 | Read paired, both mapped, this is read 1, forward strand |
| 147 | Read paired, both mapped, this is read 2, reverse strand |
| 4 | Read unmapped |
| 1796 | Read is a PCR duplicate |

```bash
# Decode any FLAG value using SAMtools
samtools flags 99
# Output: 0x63    99    PAIRED,PROPER_PAIR,MATE_REVERSE,READ1
```

### BAM Format

BAM is simply the **binary compressed version of SAM**. It is not human-readable but is far smaller and allows fast random access when indexed. Always work with BAM files in practice — SAM files are too large for real data.

```bash
# Convert SAM to BAM
samtools view -bS aligned.sam -o aligned.bam

# In practice, pipe directly from BWA-MEM2 to SAMtools to avoid creating a SAM file at all
bwa-mem2 mem ref.fa R1.fastq.gz R2.fastq.gz | samtools view -bS - -o aligned.bam
```

---

## 5.5 Post-Alignment Processing

After alignment, the BAM file needs several processing steps before it is ready for variant calling. These steps are **mandatory** for GATK Best Practices variant calling.

### Step 1 — Add Read Group Information

Read groups are metadata tags embedded in the BAM file that identify which sample, library, and sequencing run a read came from. **GATK requires read group information** — it will refuse to run without it.

Read group tags:
- `ID` — unique identifier for this read group (usually flowcell + lane)
- `SM` — sample name (the biological sample)
- `PL` — platform (ILLUMINA)
- `LB` — library name (the specific library preparation)
- `PU` — platform unit (flowcell barcode + lane)

Read groups are added during alignment using the `-R` flag:

```bash
bwa-mem2 mem \
    -R "@RG\tID:SRR6821753\tSM:SAMPLE1\tPL:ILLUMINA\tLB:lib1\tPU:SRR6821753" \
    chr20.fa \
    trimmed/sample_R1_trimmed.fastq.gz \
    trimmed/sample_R2_trimmed.fastq.gz \
    | samtools view -bS - -o aligned/sample_aligned.bam
```

### Step 2 — Sort by Coordinate

The raw output of BWA-MEM2 is in the order reads were processed — essentially random with respect to genomic position. **GATK requires coordinate-sorted BAM files.**

```bash
samtools sort -@ 8 -o aligned/sample_sorted.bam aligned/sample_aligned.bam
```

### Step 3 — Index the BAM File

An index file (`.bai`) allows tools to jump directly to any genomic region in the BAM without reading the whole file — essential for visualisation in IGV and for GATK.

```bash
samtools index aligned/sample_sorted.bam
# Creates: sample_sorted.bam.bai
```

---

## 5.6 Marking Duplicate Reads

### Why Mark Duplicates?

As discussed in Module 2, PCR duplicates are multiple sequenced copies of the same original DNA molecule, introduced during library amplification. Because they are not independent observations of the genome, they can:

- Artificially inflate the apparent evidence for a variant (if the duplicate happens to carry a sequencing error, that error appears multiple times)
- Distort allele frequencies — the ratio of reference to alternate allele reads at a variant site

**Picard MarkDuplicates** (now part of GATK4) identifies duplicate reads by finding all read pairs that have identical 5' mapping coordinates on both reads. It flags the duplicates in the BAM file (setting the duplicate FLAG bit) but does not remove them — most downstream tools then ignore flagged duplicates automatically.

### Running GATK MarkDuplicates

```bash
gatk MarkDuplicates \
    -I aligned/sample_sorted.bam \
    -O aligned/sample_markdup.bam \
    -M reports/sample_markdup_metrics.txt \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

# Re-index the marked duplicates BAM
samtools index aligned/sample_markdup.bam
```

`--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500` — on patterned flow cells (NovaSeq), optical duplicates (clusters that are physically adjacent on the flow cell and misidentified as separate) can be separated by up to 2500 pixels. Use 100 for older unpatterned flow cells (HiSeq 2500, MiSeq).

### Reading the Duplicate Metrics File

```bash
cat reports/sample_markdup_metrics.txt
```

Key fields to check:
- **PERCENT_DUPLICATION** — proportion of reads that are duplicates; > 20% for WGS warrants investigation
- **ESTIMATED_LIBRARY_SIZE** — estimated number of unique molecules in the library; useful for assessing library complexity

---

## 5.7 Base Quality Score Recalibration (BQSR)

### What is BQSR?

Base Quality Score Recalibration is a GATK preprocessing step that corrects systematic errors in the Phred quality scores assigned by the sequencer's base caller. The base caller assigns quality scores based on optical signal intensity, but these scores are not perfectly calibrated — they tend to be slightly over- or under-estimated in a systematic way that varies by:

- Cycle number (position in the read)
- Dinucleotide context (the base before the current base)
- Machine and read group

BQSR builds a model of these systematic errors by comparing the reported quality scores against the actual observed mismatch rates at known variant sites (using a database like dbSNP to exclude genuine variants from the error calculation). It then adjusts the quality scores in the BAM file to better reflect true error probabilities.

**BQSR is recommended for germline and somatic variant calling with GATK** but is not strictly required for all applications. It has the most impact on indel calling and low-frequency somatic variant detection.

### Running BQSR — Two Steps

**Step 1 — Build the recalibration model:**

```bash
# Download known variant sites (dbSNP) for BQSR
# For the workshop, use a chr20-subset of dbSNP
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi

gatk BaseRecalibrator \
    -I aligned/sample_markdup.bam \
    -R reference/chr20.fa \
    --known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -O reports/sample_recal.table
```

**Step 2 — Apply the recalibration:**

```bash
gatk ApplyBQSR \
    -I aligned/sample_markdup.bam \
    -R reference/chr20.fa \
    --bqsr-recal-file reports/sample_recal.table \
    -O aligned/sample_bqsr.bam

# Index the final BAM
samtools index aligned/sample_bqsr.bam
```

The output `sample_bqsr.bam` is your **analysis-ready BAM** — the final product of the alignment pipeline and the direct input for variant calling in Module 6.

---

## 5.8 Alignment Quality Metrics

Before proceeding to variant calling, always verify that the alignment looks as expected. Three SAMtools commands give you a comprehensive picture.

### samtools flagstat

Gives a quick count of reads in each alignment category:

```bash
samtools flagstat aligned/sample_bqsr.bam
```

Example output:
```
47382951 + 0 in total (QC-passed reads + QC-failed reads)
47382951 + 0 primary
0 + 0 secondary
0 + 0 supplementary
4821023 + 0 duplicates
42561928 + 0 primary duplicates
46891203 + 0 mapped (98.97% : N/A)
46891203 + 0 primary mapped (98.97% : N/A)
47382951 + 0 paired in sequencing
23691475 + 0 read1
23691476 + 0 read2
46012341 + 0 with itself and mate mapped
878862 + 0 singletons (1.85% : N/A)
```

**What to look for:**
- **% mapped** — should be > 95% for a clean human WGS sample; lower values suggest contamination, wrong reference, or poor library quality
- **% duplicates** — see Module 5.6; flag if > 20% for WGS
- **% singletons** — reads whose mate is unmapped; ideally < 5%

### samtools idxstats

Shows the number of reads mapped to each chromosome:

```bash
samtools idxstats aligned/sample_bqsr.bam
```

Useful for spotting unexpected mapping to non-human sequences or uneven chromosome coverage that might indicate a contaminated sample.

### samtools coverage / mosdepth

For coverage analysis — the average depth across the genome:

```bash
# Quick coverage summary with SAMtools
samtools coverage aligned/sample_bqsr.bam

# More detailed coverage with mosdepth (install separately)
conda install -c bioconda mosdepth
mosdepth --threads 8 --quantize 0:5:10:30:100: \
    reports/sample_coverage \
    aligned/sample_bqsr.bam
```

Key coverage metrics:
- **Mean depth** — average reads per base; should meet your target (e.g. 30× for WGS)
- **% bases ≥ 20×** — percentage of the genome covered at ≥ 20× depth; a standard clinical adequacy threshold
- **Coverage uniformity** — how even the coverage is across the genome; GC bias causes peaks and troughs

### GATK CollectAlignmentSummaryMetrics

The most comprehensive alignment QC for GATK pipelines:

```bash
gatk CollectAlignmentSummaryMetrics \
    -I aligned/sample_bqsr.bam \
    -R reference/chr20.fa \
    -O reports/sample_alignment_metrics.txt
```

---

## 5.9 The Full Alignment Pipeline — Summary

Putting it all together, the complete alignment pipeline from trimmed FASTQ to analysis-ready BAM is:

```
Trimmed FASTQ (from Module 4)
        ↓
1. BWA-MEM2 align + add read groups
        ↓
2. SAMtools sort (coordinate sort)
        ↓
3. SAMtools index
        ↓
4. GATK MarkDuplicates
        ↓
5. SAMtools index (re-index after markdup)
        ↓
6. GATK BaseRecalibrator (build BQSR model)
        ↓
7. GATK ApplyBQSR
        ↓
8. SAMtools index (final index)
        ↓
Analysis-ready BAM  →  Module 6: Variant Calling
```

---

## 5.10 Hands-on Exercises

### Setup

```bash
# Create directory structure for Module 5
mkdir -p ~/ngs_workshop/module5/{aligned,reference,reports,logs}
cd ~/ngs_workshop/module5

# Install tools
conda install -c bioconda bwa-mem2 samtools gatk4 mosdepth -y

# Symlink the trimmed data from Module 4
ln -s ~/ngs_workshop/module4/trimmed ./trimmed
```

---

### Exercise 5.1 — Download and Prepare the Reference

```bash
cd ~/ngs_workshop/module5/reference

# Download chr20 reference
wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
gunzip chr20.fa.gz

echo "Reference downloaded: $(wc -c < chr20.fa) bytes"

# Index with BWA-MEM2
echo "Indexing reference — this takes ~3 minutes..."
bwa-mem2 index chr20.fa

# Create FASTA index and sequence dictionary
samtools faidx chr20.fa
gatk CreateSequenceDictionary -R chr20.fa

echo "Reference ready!"
ls -lh chr20.fa*
```

---

### Exercise 5.2 — Align Reads with BWA-MEM2

```bash
cd ~/ngs_workshop/module5

# Run alignment — pipe directly to SAMtools sort to save disk space
bwa-mem2 mem \
    -t 4 \
    -R "@RG\tID:SRR6821753\tSM:SAMPLE1\tPL:ILLUMINA\tLB:lib1\tPU:SRR6821753.1" \
    reference/chr20.fa \
    trimmed/SRR6821753_R1_trimmed.fastq.gz \
    trimmed/SRR6821753_R2_trimmed.fastq.gz \
    2> logs/bwa_align.log \
    | samtools sort -@ 4 -o aligned/sample_sorted.bam -

# Index the sorted BAM
samtools index aligned/sample_sorted.bam

echo "Alignment complete"
ls -lh aligned/
```

**Check the alignment log:**
```bash
tail -5 logs/bwa_align.log
```

**Questions:**
1. How long did alignment take?
2. What does the log report about processed read pairs?

---

### Exercise 5.3 — Check Alignment Quality with flagstat

```bash
samtools flagstat aligned/sample_sorted.bam | tee reports/sample_flagstat.txt
```

Fill in the table:

```
Metric                    | Your result | Expected (good)
--------------------------|-------------|----------------
Total reads               |             | ~500,000
% mapped                  |             | > 95%
% properly paired         |             | > 90%
% duplicates (raw)        |             | varies
% singletons              |             | < 5%
```

---

### Exercise 5.4 — Mark Duplicates

```bash
gatk MarkDuplicates \
    -I aligned/sample_sorted.bam \
    -O aligned/sample_markdup.bam \
    -M reports/sample_markdup_metrics.txt \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    2> logs/markdup.log

samtools index aligned/sample_markdup.bam

# View duplicate metrics
grep -A 2 "PERCENT_DUPLICATION" reports/sample_markdup_metrics.txt
```

**Questions:**
1. What is the duplicate rate for this sample?
2. What is the estimated library size?
3. Is this duplication rate acceptable for WGS?

---

### Exercise 5.5 — Run BQSR

```bash
# Download known sites VCF for chr20
wget -q -O reference/dbsnp138_chr20.vcf.gz \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
wget -q -O reference/dbsnp138_chr20.vcf.gz.tbi \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi

# Step 1 — Build recalibration model
gatk BaseRecalibrator \
    -I aligned/sample_markdup.bam \
    -R reference/chr20.fa \
    --known-sites reference/dbsnp138_chr20.vcf.gz \
    -O reports/sample_recal.table \
    2> logs/bqsr_step1.log

# Step 2 — Apply BQSR
gatk ApplyBQSR \
    -I aligned/sample_markdup.bam \
    -R reference/chr20.fa \
    --bqsr-recal-file reports/sample_recal.table \
    -O aligned/sample_bqsr.bam \
    2> logs/bqsr_step2.log

samtools index aligned/sample_bqsr.bam

echo "Analysis-ready BAM created: aligned/sample_bqsr.bam"
ls -lh aligned/
```

---

### Exercise 5.6 — Final BAM Quality Check

```bash
# Run all QC metrics on the final BAM
echo "=== flagstat ===" && samtools flagstat aligned/sample_bqsr.bam
echo ""
echo "=== coverage summary ===" && samtools coverage aligned/sample_bqsr.bam | head -5
echo ""
echo "=== idxstats ===" && samtools idxstats aligned/sample_bqsr.bam
```

---

### Exercise 5.7 — Write the Full Pipeline as a Single Script

```bash
cat > ~/ngs_workshop/module5/run_alignment.sh << 'EOF'
#!/bin/bash
# Full alignment pipeline: trimmed FASTQ → analysis-ready BAM
# Usage: bash run_alignment.sh <sample_name> <R1_trimmed> <R2_trimmed>

set -euo pipefail

SAMPLE=$1
R1=$2
R2=$3

REF=reference/chr20.fa
KNOWN_SITES=reference/dbsnp138_chr20.vcf.gz
THREADS=4

mkdir -p aligned reports logs

echo "[$(date)] Starting alignment pipeline for: ${SAMPLE}"

# Step 1 — Align
echo "[$(date)] Step 1/5: Aligning with BWA-MEM2..."
bwa-mem2 mem -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:${SAMPLE}.1" \
    ${REF} ${R1} ${R2} 2> logs/${SAMPLE}_bwa.log \
    | samtools sort -@ ${THREADS} -o aligned/${SAMPLE}_sorted.bam -
samtools index aligned/${SAMPLE}_sorted.bam

# Step 2 — Mark duplicates
echo "[$(date)] Step 2/5: Marking duplicates..."
gatk MarkDuplicates \
    -I aligned/${SAMPLE}_sorted.bam \
    -O aligned/${SAMPLE}_markdup.bam \
    -M reports/${SAMPLE}_markdup_metrics.txt \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    2> logs/${SAMPLE}_markdup.log
samtools index aligned/${SAMPLE}_markdup.bam

# Step 3 — BQSR
echo "[$(date)] Step 3/5: Running BQSR..."
gatk BaseRecalibrator \
    -I aligned/${SAMPLE}_markdup.bam \
    -R ${REF} \
    --known-sites ${KNOWN_SITES} \
    -O reports/${SAMPLE}_recal.table \
    2> logs/${SAMPLE}_bqsr1.log

gatk ApplyBQSR \
    -I aligned/${SAMPLE}_markdup.bam \
    -R ${REF} \
    --bqsr-recal-file reports/${SAMPLE}_recal.table \
    -O aligned/${SAMPLE}_bqsr.bam \
    2> logs/${SAMPLE}_bqsr2.log
samtools index aligned/${SAMPLE}_bqsr.bam

# Step 4 — QC metrics
echo "[$(date)] Step 4/5: Collecting QC metrics..."
samtools flagstat aligned/${SAMPLE}_bqsr.bam > reports/${SAMPLE}_flagstat.txt
samtools coverage aligned/${SAMPLE}_bqsr.bam > reports/${SAMPLE}_coverage.txt

# Step 5 — Cleanup intermediate files
echo "[$(date)] Step 5/5: Cleaning up intermediate files..."
rm -f aligned/${SAMPLE}_sorted.bam aligned/${SAMPLE}_sorted.bam.bai
rm -f aligned/${SAMPLE}_markdup.bam aligned/${SAMPLE}_markdup.bam.bai

echo "[$(date)] Pipeline complete!"
echo "Analysis-ready BAM: aligned/${SAMPLE}_bqsr.bam"
echo "QC report: reports/${SAMPLE}_flagstat.txt"
EOF

chmod +x ~/ngs_workshop/module5/run_alignment.sh

# Run the pipeline
bash ~/ngs_workshop/module5/run_alignment.sh \
    SRR6821753 \
    trimmed/SRR6821753_R1_trimmed.fastq.gz \
    trimmed/SRR6821753_R2_trimmed.fastq.gz
```

---

> 💡 **Key takeaway for this module:** The alignment pipeline produces far more than just a BAM file — each step (read groups, sorting, duplicate marking, BQSR) exists for a specific reason tied to the requirements of downstream variant calling. Skipping any step will either cause GATK to fail outright or silently degrade the quality of your variant calls. The analysis-ready BAM you produce here — `sample_bqsr.bam` — is the direct input for Module 6.

---

**Previous:** [Module 4 — Data Preprocessing and Cleaning](./module4_preprocessing.md)  
**Next:** [Module 6 — Variant Detection and Analysis](./module6_variant_calling.md)
