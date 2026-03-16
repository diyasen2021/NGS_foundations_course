# Module 7: Visualization and Interpretation

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 2 | Module 7 of 4 | 🧪 Hands-on

> **Context:** You now have an annotated VCF file from Module 6 and an analysis-ready BAM from Module 5. This module covers how to visually inspect your alignments and variants, interpret what you see, and communicate findings clearly.

---

## Table of Contents

- [7.1 Why Visualisation Matters](#71-why-visualisation-matters)
- [7.2 IGV — Integrative Genomics Viewer](#72-igv--integrative-genomics-viewer)
- [7.3 Loading Data into IGV](#73-loading-data-into-igv)
- [7.4 Interpreting Alignments in IGV](#74-interpreting-alignments-in-igv)
- [7.5 Visualising Variants in IGV](#75-visualising-variants-in-igv)
- [7.6 Recognising Artefacts vs Real Variants in IGV](#76-recognising-artefacts-vs-real-variants-in-igv)
- [7.7 UCSC Genome Browser](#77-ucsc-genome-browser)
- [7.8 Coverage Plots and Summary Visualisations in R](#78-coverage-plots-and-summary-visualisations-in-r)
- [7.9 Visualising Variant Annotations](#79-visualising-variant-annotations)
- [7.10 Hands-on Exercises](#710-hands-on-exercises)

---

## 7.1 Why Visualisation Matters

After running a bioinformatics pipeline, it is tempting to trust the numbers — mapping rate 98.7%, 4,321 PASS variants, Ti/Tv ratio 2.1 — and move straight to biological interpretation. This is a mistake. Every experienced genomics analyst knows that pipelines produce artefacts, and that manual visualisation of key findings is a non-negotiable step before drawing any conclusions.

Visualisation allows you to:

- **Manually confirm variants** — see with your own eyes whether the reads supporting a variant look genuine or like an alignment artefact
- **Identify regional problems** — spot regions of low coverage, repetitive sequence, or systematic misalignment that automated tools may not flag
- **Understand your data** — develop an intuition for what good and bad data looks like, which makes you a better analyst
- **Communicate findings** — generate figures that explain your results to collaborators, supervisors, and in publications

The two most important visualisation tools in NGS are **IGV** (for looking at reads and variants at specific genomic loci) and the **UCSC Genome Browser** (for contextualising your findings within genome annotation tracks). Both are covered in detail in this module.

---

## 7.2 IGV — Integrative Genomics Viewer

### What is IGV?

The Integrative Genomics Viewer (IGV) was developed at the Broad Institute and is the most widely used desktop tool for visualising NGS alignment data. It allows you to load BAM files, VCF files, BED files, and annotation tracks, and navigate the genome interactively to inspect alignments at any region.

IGV is free, runs on Windows, Mac, and Linux, and can handle very large files efficiently by only loading the data in the current view window rather than the entire file at once — which is why BAM files must be indexed (`.bai`) before loading.

### Installation

```bash
# Download IGV from the Broad Institute website
# https://igv.org/doc/desktop/#DownloadPage/

# On Linux — download and extract
wget https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_Linux_2.17.4_WithJava.zip
unzip IGV_Linux_2.17.4_WithJava.zip
cd IGV_Linux_2.17.4
./igv.sh

# On macOS — download the .app bundle from https://igv.org
# On Windows — download the Windows installer from https://igv.org

# Alternative — IGV web app (no installation required)
# Open in browser: https://igv.org/app
```

> 💡 **For GitHub Codespaces:** Since Codespaces is a headless Linux environment, the desktop IGV application cannot run directly. Use the **IGV web app** at https://igv.org/app instead — it runs entirely in the browser and accepts the same BAM and VCF files. Alternatively, generate IGV screenshots using igvtools (covered in Section 7.10).

### IGV Interface Overview

```
┌─────────────────────────────────────────────────────┐
│  Toolbar: genome selector, navigation, zoom          │
├─────────────────────────────────────────────────────┤
│  Chromosome ideogram — click to jump to region       │
├─────────────────────────────────────────────────────┤
│  Gene annotation track (RefSeq genes)                │
├─────────────────────────────────────────────────────┤
│  Coverage track — depth histogram                    │
├─────────────────────────────────────────────────────┤
│  Read alignment track — individual reads as arrows   │
├─────────────────────────────────────────────────────┤
│  VCF track — variant calls as coloured marks         │
└─────────────────────────────────────────────────────┘
```

**Key navigation shortcuts:**

| Action | How |
|---|---|
| Go to a gene or position | Type in the search box: `chr20:1,234,567` or `BRCA1` |
| Zoom in | `+` key or scroll wheel up |
| Zoom out | `-` key or scroll wheel down |
| Pan left/right | Click and drag, or arrow keys |
| Jump to next variant | Click VCF track; use arrow keys |
| Increase track height | Drag bottom edge of track downward |
| Sort reads | Right-click on read track → Sort reads by |
| Colour reads | Right-click on read track → Color reads by |

---

## 7.3 Loading Data into IGV

### What Files to Load

For this workshop, you will load three types of files into IGV:

| File | What it shows |
|---|---|
| `sample_bqsr.bam` + `sample_bqsr.bam.bai` | Individual read alignments |
| `sample_PASS.vcf.gz` + `.tbi` | Variant call positions |
| RefSeq gene annotation | Gene models (built into IGV) |

### Loading Files — Step by Step

**1. Set the reference genome**
- Open IGV
- Click the genome dropdown (top left) → select **Human (hg38)**
- This loads the GRCh38 gene annotation automatically

**2. Load the BAM file**
- File → Load from File → select `sample_bqsr.bam`
- IGV automatically finds the `.bai` index file in the same directory
- Two tracks appear: a **coverage track** (grey histogram) and a **read track** (coloured arrows)

**3. Load the VCF file**
- File → Load from File → select `sample_PASS.vcf.gz`
- A new track appears showing variant positions as coloured marks

**4. Navigate to a region of interest**
- Type a position in the search box: `chr20:31,022,000-31,025,000`
- Or type a gene name: `RUNX1T1` (a gene on chr20)

### Loading Files for IGV Web App

```
1. Open https://igv.org/app
2. Genome → hg38
3. Tracks → Local File → select sample_bqsr.bam
   (also select sample_bqsr.bam.bai at the same time)
4. Tracks → Local File → select sample_PASS.vcf.gz
   (also select sample_PASS.vcf.gz.tbi)
```

---

## 7.4 Interpreting Alignments in IGV

### The Coverage Track

The grey histogram at the top of the read track shows the **read depth** at each position — the number of reads covering that base. A consistent, even coverage histogram is a good sign. Look for:

- **Coverage gaps** (depth = 0) — regions with no reads; variants here cannot be called
- **Coverage spikes** (very high depth) — may indicate repetitive regions, PCR amplification bias, or mapping of reads from multiple genomic locations to one spot
- **Coloured bars in the coverage track** — appear when a base differs from the reference in > 20% of reads at that position; this is how IGV flags potential variant sites visually

### Reading Individual Reads

Each read is shown as a grey arrow pointing in the direction of sequencing (right-pointing = forward strand, left-pointing = reverse strand). Read pairs are connected by a thin line.

**Read colours and what they mean:**

| Colour | Meaning |
|---|---|
| Grey | Normal read — maps to expected location with expected insert size |
| White | Low mapping quality (MAPQ < 20) |
| Red | Insert size larger than expected — may indicate deletion |
| Blue | Insert size smaller than expected — may indicate insertion |
| Coloured base | A base that differs from the reference |

**Base colours at variant positions:**

| Colour | Base |
|---|---|
| Green | A |
| Blue | C |
| Orange/Brown | G |
| Red | T |
| Grey | Deletion |
| Purple | Insertion |

### Coverage Depth and Allele Balance

At a heterozygous variant site (GT = 0/1), you expect approximately half the reads to carry the reference allele and half to carry the alternate allele. In IGV, you should see roughly equal numbers of grey (reference) and coloured (alternate) bases stacked at that position.

A healthy heterozygous SNP looks like:

```
Coverage:  ████████████████████████  (uniform depth ~30×)
           AAAAAAAAGGGGGGGGGAAAAAA
           AAAAAAAAGGGGGGGGGAAAAAA   ← ~50% reads show G (alternate)
           AAAAAAAAGGGGGGGGGAAAAAA
           AAAAAAAAAAAAAAAAAAAAAAA   ← ~50% reads show A (reference)
           AAAAAAAAAAAAAAAAAAAAAAA
```

---

## 7.5 Visualising Variants in IGV

### SNPs

A true heterozygous SNP in IGV shows:
- Coloured bases at the variant position in approximately 50% of reads
- Both strands (forward and reverse arrows) carry the variant equally
- The variant appears consistently across the full width of reads — not just at read ends
- No unusual insert sizes or soft-clipping around the site

A true homozygous SNP shows:
- 100% of reads carry the alternate allele — no grey (reference) bases visible at that position

### Insertions

Insertions appear as **purple I** markers in IGV. Clicking on the marker shows the inserted sequence. The reads flanking the insertion should show clean alignment on both sides.

```
Reference:  ATCGAT----GCTAGC
Reads:      ATCGATTTTTGCTAGC  ← 4 bp insertion (shown as purple I)
            ATCGATTTTTGCTAGC
            ATCGAT----GCTAGC  ← reads without insertion (reference allele)
```

### Deletions

Deletions appear as **grey horizontal lines** (gaps) in the reads. The reads show continuous alignment on both sides of the deletion, with the deleted bases simply absent.

```
Reference:  ATCGATTTTTGCTAGC
Reads:      ATCGAT----GCTAGC  ← 4 bp deletion (shown as gap/line)
            ATCGAT----GCTAGC
            ATCGATTTTTGCTAGC  ← reads with reference allele
```

### What to Check for Every Important Variant

Before reporting or acting on any variant, inspect it in IGV and confirm:

```
□ Does the variant appear on both forward and reverse strand reads?
  → If only on one strand — possible strand bias artefact

□ Does the variant appear throughout the read, not just at ends?
  → Variants only at read ends suggest low-quality base calls

□ Is the depth sufficient? (ideally ≥ 20× for germline calling)
  → Low depth means the call is uncertain

□ Is the allele balance reasonable?
  → Het variant: expect ~40-60% ALT reads
  → Hom variant: expect ≥ 90% ALT reads

□ Are the surrounding reads well-aligned (no excessive soft-clipping)?
  → Excessive soft-clips nearby suggest a difficult/repetitive region

□ Is the mapping quality high (reads shown in normal grey, not white)?
  → White reads have MAPQ < 20 and should not be trusted
```

---

## 7.6 Recognising Artefacts vs Real Variants in IGV

Learning to distinguish genuine variants from artefacts is one of the most valuable skills in clinical and research genomics. Here are the most common artefact patterns:

### Strand Bias Artefact

A variant supported only by reads on one strand (all forward or all reverse) is suspicious. True variants arise from the original DNA and should be sequenced on both strands during library preparation.

**In IGV:** Right-click read track → Color reads by → First-of-pair strand. If all ALT reads are one colour (all forward or all reverse), this is a strand bias artefact.

```
Forward reads:  ATCG[G]GCTA  ← ALT allele only on forward strand
Forward reads:  ATCG[G]GCTA
Reverse reads:  ATCG[A]GCTA  ← only reference on reverse strand
Reverse reads:  ATCG[A]GCTA
→ SUSPICIOUS — likely artefact
```

### End-of-Read Artefact

Variants that only appear in the last 5–10 bases of reads are often low-quality base calls rather than real variants. Base quality declines at read ends.

**In IGV:** The coloured bases appear only at the very tips of reads, never in the middle.

### Misalignment in Repetitive Regions

In repetitive regions, reads from elsewhere in the genome can be incorrectly mapped, creating phantom variant signals. These regions typically show:
- Very high or irregular coverage
- Reads with low MAPQ (shown as white in IGV)
- Variants that do not follow expected allele balance

**In IGV:** Right-click → Color reads by → Mapping quality. White reads = MAPQ 0 = unreliable. If most reads at a variant site are white, the call should not be trusted.

### FFPE Artefact (C→T / G→A)

In FFPE-derived samples, formalin fixation causes cytosine deamination — C→T changes (and G→A on the opposite strand). These appear as real SNPs but are chemical artefacts.

**In IGV:** C→T variants in FFPE samples should show strong strand bias — deamination preferentially affects one strand. Check with the strand bias filter (FS or SOR in the VCF).

### PCR Slippage (Stutter) at Indels

In repetitive sequence (e.g. a run of 10 adenines), PCR polymerase can slip during amplification, creating indels of 1–2 bp that are artefactual. These appear as low-frequency indels clustered in simple repeats.

**In IGV:** The indel occurs in the middle of a homopolymer or short tandem repeat. Check the surrounding sequence — if the region is repetitive, treat the indel call with scepticism.

---

## 7.7 UCSC Genome Browser

### What is the UCSC Genome Browser?

The UCSC Genome Browser (https://genome.ucsc.edu) is a web-based tool for exploring the genome in the context of hundreds of annotation tracks — gene models, conservation scores, regulatory elements, repeat elements, population variants, and more. Unlike IGV, which focuses on your own data, the UCSC browser is primarily used for contextualising variants within existing biological knowledge.

### Key Uses in Variant Interpretation

**1. Looking up a gene or variant**  
Type a gene name, genomic coordinates, or rsID into the search box. The browser shows the gene structure, known transcripts, nearby genes, and any overlapping annotations.

**2. Checking if a region is repetitive**  
Enable the **RepeatMasker** track to see if your variant falls within a SINE, LINE, or satellite repeat. Variants in repetitive regions should be treated with caution.

**3. Checking conservation**  
The **PhyloP** and **PhastCons** tracks show evolutionary conservation across vertebrates. A variant in a highly conserved position is more likely to be functionally important.

**4. Checking regulatory context**  
The **ENCODE** regulation tracks show DNase hypersensitive sites (open chromatin), transcription factor binding sites, and histone modification signals. A variant in an active regulatory element may affect gene expression even if it is in a non-coding region.

**5. Checking population variant databases**  
The **gnomAD** and **dbSNP** tracks show known common variants. If your variant of interest overlaps a known common variant (AF > 1%), it is likely a benign polymorphism.

### Navigating to Your Variant

```
1. Go to https://genome.ucsc.edu
2. Click "Genome Browser" → select Human GRCh38/hg38
3. In the search box type your variant coordinates:
   chr20:31,022,459  or  rs123456789
4. Click "Go"
5. The browser shows the region with default annotation tracks
```

### Adding Custom Tracks

You can upload your own VCF file directly to the UCSC browser as a custom track:

```
1. In the browser, click "My Data" → "Custom Tracks"
2. Click "Choose File" → select your sample_PASS.vcf.gz
3. Click "Submit"
4. Your variants appear as a new track in the browser
```

---

## 7.8 Coverage Plots and Summary Visualisations in R

While IGV and the UCSC browser are excellent for inspecting individual loci, summary visualisations help you understand the overall quality and characteristics of your dataset. R and Python are the standard tools for generating these plots.

### Setting Up R

```bash
# Install R and required packages in Codespaces
conda install -c conda-forge r-base r-ggplot2 r-dplyr r-tidyr -y

# Open R
R
```

### Plot 1 — Coverage Distribution

```r
library(ggplot2)
library(dplyr)

# Read coverage data generated by samtools coverage (Module 5)
cov <- read.table("reports/sample_coverage.txt",
                  header = TRUE, sep = "\t", comment.char = "#")

# Coverage histogram
ggplot(cov, aes(x = meandepth)) +
  geom_histogram(bins = 50, fill = "#0D7377", colour = "white") +
  geom_vline(xintercept = 30, linetype = "dashed", colour = "red", linewidth = 0.8) +
  labs(
    title = "Coverage Distribution per Chromosome",
    x = "Mean Depth (×)",
    y = "Number of Regions",
    caption = "Red dashed line = 30× target coverage"
  ) +
  theme_minimal(base_size = 13)

ggsave("reports/coverage_distribution.png", width = 8, height = 5, dpi = 150)
```

### Plot 2 — Variant Quality Score Distribution

```r
# Read the variants table generated by GATK VariantsToTable
variants <- read.table("variants/sample_variants_table.tsv",
                       header = TRUE, sep = "\t", na.strings = ".")

# Quality score distribution
ggplot(variants, aes(x = QUAL)) +
  geom_histogram(bins = 60, fill = "#14A085", colour = "white") +
  geom_vline(xintercept = 30, linetype = "dashed", colour = "red", linewidth = 0.8) +
  scale_x_log10() +
  labs(
    title = "Variant Quality Score Distribution",
    x = "QUAL score (log scale)",
    y = "Number of Variants"
  ) +
  theme_minimal(base_size = 13)

ggsave("reports/variant_quality_distribution.png", width = 8, height = 5, dpi = 150)
```

### Plot 3 — Allele Depth (AD) — Checking Allele Balance

```r
# Parse AD field — it is stored as "ref_depth,alt_depth"
variants_ad <- variants %>%
  filter(!is.na(SAMPLE1.AD)) %>%
  tidyr::separate(SAMPLE1.AD, into = c("ref_depth", "alt_depth"),
                  sep = ",", convert = TRUE) %>%
  mutate(
    total_depth = ref_depth + alt_depth,
    vaf = alt_depth / total_depth
  ) %>%
  filter(total_depth > 0)

# VAF distribution
ggplot(variants_ad, aes(x = vaf)) +
  geom_histogram(bins = 50, fill = "#0A2342", colour = "white") +
  geom_vline(xintercept = c(0.5, 1.0), linetype = "dashed",
             colour = c("#F4A261", "#E63946"), linewidth = 0.8) +
  labs(
    title = "Variant Allele Frequency Distribution",
    x = "Variant Allele Frequency (VAF)",
    y = "Number of Variants",
    caption = "Orange = 0.5 (expected het); Red = 1.0 (expected hom)"
  ) +
  theme_minimal(base_size = 13)

ggsave("reports/vaf_distribution.png", width = 8, height = 5, dpi = 150)
```

### Plot 4 — Variant Type Breakdown

```r
# Count variant types from SnpEff annotation
# Parse the ANN field from the snpEff annotated VCF
library(vcfR)

# Read the annotated VCF
vcf <- read.vcfR("variants/sample_snpeff.vcf.gz", verbose = FALSE)

# Extract ANN field and parse effect
ann <- extract.info(vcf, "ANN")
effects <- sapply(strsplit(ann, "\\|"), function(x) x[2])

effect_counts <- as.data.frame(table(effects)) %>%
  arrange(desc(Freq)) %>%
  head(15)

ggplot(effect_counts, aes(x = reorder(effects, Freq), y = Freq)) +
  geom_col(fill = "#0D7377") +
  coord_flip() +
  labs(
    title = "Variant Effects — Top 15 Categories",
    x = "SnpEff Effect",
    y = "Number of Variants"
  ) +
  theme_minimal(base_size = 12)

ggsave("reports/variant_effects.png", width = 9, height = 6, dpi = 150)
```

### Plot 5 — Ti/Tv Ratio per Chromosome

```r
# Read the raw variants table
snps <- variants %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1)  # SNPs only

# Classify as transition or transversion
transitions <- c("AG", "GA", "CT", "TC")
snps <- snps %>%
  mutate(
    change = paste0(REF, ALT),
    type = ifelse(change %in% transitions, "Transition", "Transversion")
  )

titv <- snps %>%
  group_by(CHROM, type) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = type, values_from = count, values_fill = 0) %>%
  mutate(TiTv = Transition / Transversion)

ggplot(titv, aes(x = CHROM, y = TiTv)) +
  geom_col(fill = "#14A085") +
  geom_hline(yintercept = 2.1, linetype = "dashed", colour = "red", linewidth = 0.8) +
  labs(
    title = "Ti/Tv Ratio per Chromosome",
    x = "Chromosome",
    y = "Ti/Tv Ratio",
    caption = "Red line = expected WGS Ti/Tv ~2.1"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("reports/titv_per_chromosome.png", width = 10, height = 5, dpi = 150)
```

---

## 7.9 Visualising Variant Annotations

### Filtering and Exploring the Annotated VCF in R

```r
# Read the variants table
variants <- read.table("variants/sample_variants_table.tsv",
                       header = TRUE, sep = "\t", na.strings = ".")

# Summary of genotype counts
table(variants$SAMPLE1.GT)

# How many variants have GQ >= 20?
sum(variants$SAMPLE1.GQ >= 20, na.rm = TRUE)

# Distribution of read depth at variant sites
summary(variants$SAMPLE1.DP)
```

### Using bcftools to Filter and Summarise

```bash
# Count variants by chromosome
bcftools stats -r chr20 variants/sample_PASS.vcf.gz | grep "^PSC"

# Extract all missense variants to a separate file
bcftools view variants/sample_snpeff.vcf.gz | \
    grep -E "^#|missense_variant" | \
    bgzip > variants/sample_missense.vcf.gz
tabix -p vcf variants/sample_missense.vcf.gz

# Count high-impact variants (stop gained, frameshift)
bcftools view variants/sample_snpeff.vcf.gz | \
    grep -cE "stop_gained|frameshift_variant"

# View variants in a specific gene
bcftools view variants/sample_snpeff.vcf.gz | \
    grep -E "^#|GENE_NAME_HERE"
```

### Summary Report Table in R

```r
# Generate a clean summary table of high-priority variants
priority_variants <- variants %>%
  filter(
    SAMPLE1.DP >= 20,     # minimum depth
    SAMPLE1.GQ >= 20,     # minimum genotype quality
    QUAL >= 50            # minimum variant quality
  ) %>%
  select(CHROM, POS, REF, ALT, QUAL, SAMPLE1.GT, SAMPLE1.AD, SAMPLE1.DP, SAMPLE1.GQ) %>%
  arrange(desc(QUAL))

# Print top 20
head(priority_variants, 20)

# Write to CSV for sharing
write.csv(priority_variants, "reports/priority_variants.csv",
          row.names = FALSE, quote = FALSE)
```

---

## 7.10 Hands-on Exercises

### Setup

```bash
mkdir -p ~/ngs_workshop/module7/{igv_screenshots,reports,plots}
cd ~/ngs_workshop/module7

# Install R packages
conda install -c conda-forge r-base r-ggplot2 r-dplyr r-tidyr r-vcfr -y

# Symlink outputs from previous modules
ln -s ~/ngs_workshop/module5/aligned ./aligned
ln -s ~/ngs_workshop/module6/variants ./variants
ln -s ~/ngs_workshop/module5/reports ./alignment_reports
ln -s ~/ngs_workshop/module6/reports ./variant_reports
```

---

### Exercise 7.1 — Generate IGV Screenshots with igvtools (Headless)

Since Codespaces cannot display a desktop GUI, we use `igvtools` to generate batch screenshots from the command line.

```bash
# Install igvtools
conda install -c bioconda igvtools -y

# Create a batch script for IGV
cat > igv_batch.txt << 'EOF'
new
genome hg38
load aligned/sample_bqsr.bam
load variants/sample_PASS.vcf.gz
snapshotDirectory igv_screenshots/

# Screenshot 1 — zoom out view of chr20
goto chr20
collapse
snapshot chr20_overview.png

# Screenshot 2 — navigate to a specific variant
# Replace with an actual variant position from your VCF
goto chr20:31,022,400-31,022,600
expand
sort position
snapshot variant_locus_1.png

# Screenshot 3 — another variant locus
goto chr20:44,619,500-44,619,700
sort base
snapshot variant_locus_2.png

exit
EOF

# Run igvtools in batch mode
igvtools index aligned/sample_bqsr.bam  # ensure BAM is indexed
igv.sh -b igv_batch.txt

ls -lh igv_screenshots/
```

---

### Exercise 7.2 — Manually Inspect Variants in the IGV Web App

Open https://igv.org/app in your browser and complete the following:

```
Step 1 — Set genome to hg38

Step 2 — Load your BAM file
  Tracks → Local File → aligned/sample_bqsr.bam
  (select both .bam and .bam.bai)

Step 3 — Load your VCF file
  Tracks → Local File → variants/sample_PASS.vcf.gz
  (select both .vcf.gz and .vcf.gz.tbi)

Step 4 — Navigate to your first variant
  Look at your variants table (variants/sample_variants_table.tsv)
  Pick any variant with GQ >= 30 and DP >= 20
  Type its coordinates into the search box:
  e.g.  chr20:31,022,459

Step 5 — For each variant, fill in this checklist:
```

**Variant inspection checklist:**

```
Variant: chr20:_______ REF:___ ALT:___

□ What is the read depth at this site? ___
□ What is the approximate VAF (% ALT reads)? ___
□ Is the variant on both strands? Yes / No
□ Does it appear throughout reads (not just at ends)? Yes / No
□ Are reads well-aligned (no excessive soft-clipping)? Yes / No
□ Are reads high MAPQ (grey, not white)? Yes / No
□ What is your overall assessment? Likely real / Suspicious artefact
```

---

### Exercise 7.3 — Coverage and Variant Quality Plots in R

```bash
cd ~/ngs_workshop/module7
R --no-save << 'EOF'

library(ggplot2)
library(dplyr)

# --- Plot 1: Variant depth distribution ---
variants <- read.table("variants/sample_variants_table.tsv",
                       header=TRUE, sep="\t", na.strings=".")

p1 <- ggplot(variants, aes(x=SAMPLE1.DP)) +
  geom_histogram(bins=40, fill="#0D7377", colour="white") +
  geom_vline(xintercept=20, linetype="dashed", colour="red") +
  labs(title="Read Depth at Variant Sites",
       x="Read Depth (DP)", y="Number of Variants",
       caption="Red line = 20× minimum threshold") +
  theme_minimal(base_size=13)
ggsave("plots/variant_depth.png", p1, width=8, height=5, dpi=150)
cat("Saved: plots/variant_depth.png\n")

# --- Plot 2: Genotype quality distribution ---
p2 <- ggplot(variants, aes(x=SAMPLE1.GQ)) +
  geom_histogram(bins=40, fill="#14A085", colour="white") +
  geom_vline(xintercept=20, linetype="dashed", colour="red") +
  labs(title="Genotype Quality Distribution",
       x="Genotype Quality (GQ)", y="Number of Variants",
       caption="Red line = GQ 20 minimum threshold") +
  theme_minimal(base_size=13)
ggsave("plots/genotype_quality.png", p2, width=8, height=5, dpi=150)
cat("Saved: plots/genotype_quality.png\n")

# --- Plot 3: VAF distribution ---
variants_ad <- variants %>%
  filter(!is.na(SAMPLE1.AD)) %>%
  tidyr::separate(SAMPLE1.AD, into=c("ref_depth","alt_depth"),
                  sep=",", convert=TRUE) %>%
  mutate(total=ref_depth+alt_depth, vaf=alt_depth/total) %>%
  filter(total > 0)

p3 <- ggplot(variants_ad, aes(x=vaf)) +
  geom_histogram(bins=50, fill="#0A2342", colour="white") +
  geom_vline(xintercept=c(0.5,1.0), linetype="dashed",
             colour=c("#F4A261","#E63946")) +
  labs(title="Variant Allele Frequency Distribution",
       x="VAF", y="Number of Variants") +
  theme_minimal(base_size=13)
ggsave("plots/vaf_distribution.png", p3, width=8, height=5, dpi=150)
cat("Saved: plots/vaf_distribution.png\n")

cat("All plots saved to plots/\n")
EOF
```

---

### Exercise 7.4 — Export a Priority Variant Table

```bash
R --no-save << 'EOF'
variants <- read.table("variants/sample_variants_table.tsv",
                       header=TRUE, sep="\t", na.strings=".")

# Apply quality filters and select key columns
priority <- variants[
  !is.na(variants$SAMPLE1.DP) &
  variants$SAMPLE1.DP >= 20 &
  !is.na(variants$SAMPLE1.GQ) &
  variants$SAMPLE1.GQ >= 20 &
  variants$QUAL >= 50,
]

priority <- priority[order(-priority$QUAL), ]
priority <- priority[, c("CHROM","POS","REF","ALT","QUAL",
                          "SAMPLE1.GT","SAMPLE1.AD",
                          "SAMPLE1.DP","SAMPLE1.GQ")]

cat("Total priority variants:", nrow(priority), "\n")
cat("Heterozygous:", sum(priority$SAMPLE1.GT == "0/1", na.rm=TRUE), "\n")
cat("Homozygous ALT:", sum(priority$SAMPLE1.GT == "1/1", na.rm=TRUE), "\n")

write.csv(priority, "reports/priority_variants.csv",
          row.names=FALSE, quote=FALSE)
cat("Saved: reports/priority_variants.csv\n")
EOF
```

---

### Exercise 7.5 — Look Up a Variant in the UCSC Genome Browser

Pick one of your top-ranked variants from `reports/priority_variants.csv` and investigate it using the UCSC Genome Browser.

```
1. Go to https://genome.ucsc.edu
2. Select Genomes → Human GRCh38/hg38
3. Enter your variant coordinates in the search box
4. Answer the following questions:

   □ Which gene (if any) does this variant fall in?
   □ Which transcript/exon?
   □ Is the region covered by RepeatMasker repeats?
   □ What is the PhyloP conservation score at this position?
     (Higher = more conserved across vertebrates)
   □ Does any ENCODE track show this region is regulatory?
   □ Is this position covered by a dbSNP entry?
     If yes — what is the rsID and population frequency?
```

---

### Exercise 7.6 — Build a Complete Visualisation Report

Run all plots in one script and generate a summary:

```bash
cat > ~/ngs_workshop/module7/run_visualisation.sh << 'EOF'
#!/bin/bash
set -euo pipefail

echo "[$(date)] Generating visualisation report..."

mkdir -p plots reports

R --no-save --quiet << 'REOF'
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

variants <- read.table("variants/sample_variants_table.tsv",
                       header=TRUE, sep="\t", na.strings=".")

# 1. Depth distribution
ggsave("plots/depth.png",
  ggplot(variants, aes(x=SAMPLE1.DP)) +
    geom_histogram(bins=40, fill="#0D7377", colour="white") +
    geom_vline(xintercept=20, linetype="dashed", colour="red") +
    labs(title="Read Depth at Variant Sites", x="Depth", y="Count") +
    theme_minimal(),
  width=7, height=4, dpi=150)

# 2. GQ distribution
ggsave("plots/genotype_quality.png",
  ggplot(variants, aes(x=SAMPLE1.GQ)) +
    geom_histogram(bins=40, fill="#14A085", colour="white") +
    geom_vline(xintercept=20, linetype="dashed", colour="red") +
    labs(title="Genotype Quality", x="GQ", y="Count") +
    theme_minimal(),
  width=7, height=4, dpi=150)

# 3. VAF
variants_ad <- variants %>%
  filter(!is.na(SAMPLE1.AD)) %>%
  separate(SAMPLE1.AD, into=c("ref","alt"), sep=",", convert=TRUE) %>%
  mutate(vaf=alt/(ref+alt)) %>%
  filter((ref+alt) > 0)

ggsave("plots/vaf.png",
  ggplot(variants_ad, aes(x=vaf)) +
    geom_histogram(bins=50, fill="#0A2342", colour="white") +
    geom_vline(xintercept=c(0.5,1.0), linetype="dashed",
               colour=c("#F4A261","#E63946")) +
    labs(title="Variant Allele Frequency", x="VAF", y="Count") +
    theme_minimal(),
  width=7, height=4, dpi=150)

cat("Plots complete\n")
REOF

echo "[$(date)] All plots saved to plots/"
echo "[$(date)] Visualisation report complete"
EOF

chmod +x ~/ngs_workshop/module7/run_visualisation.sh
bash ~/ngs_workshop/module7/run_visualisation.sh
```

---

> 💡 **Key takeaway for this module:** Visualisation is not an optional extra at the end of a pipeline — it is an integral part of the analysis. Every important variant should be visually inspected in IGV before it is reported. The UCSC Genome Browser contextualises your variants within decades of accumulated genome annotation. Summary plots in R give you and your collaborators an immediate, honest picture of data quality and variant characteristics. Together these tools allow you to tell a complete, trustworthy story about your sequencing data.

---

**Previous:** [Module 6 — Variant Detection and Analysis](./module6_variant_calling.md)  
**Next:** [Module 8 — Advanced NGS Pipelines and Workflow Management](./module8_pipelines.md)
