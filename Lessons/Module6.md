# Module 6: Variant Detection and Analysis

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 2 | Module 6 of 4 | 🧪 Hands-on

> **Context:** This module takes the analysis-ready BAM file produced in Module 5 and calls variants using the GATK Best Practices pipeline. The output is an annotated VCF file containing high-confidence SNPs and indels.

---

## Table of Contents

- [6.1 What is a Variant?](#61-what-is-a-variant)
- [6.2 The VCF File Format](#62-the-vcf-file-format)
- [6.3 Germline vs Somatic Variant Calling](#63-germline-vs-somatic-variant-calling)
- [6.4 The GATK HaplotypeCaller — How it Works](#64-the-gatk-haplotypecaller--how-it-works)
- [6.5 Variant Filtering — VQSR and Hard Filters](#65-variant-filtering--vqsr-and-hard-filters)
- [6.6 Variant Annotation](#66-variant-annotation)
- [6.7 Interpreting Annotated Variants](#67-interpreting-annotated-variants)
- [6.8 Hands-on Exercises](#68-hands-on-exercises)
- [6.9 Common Problems and Troubleshooting](#69-common-problems-and-troubleshooting)
- [6.10 The Full Variant Calling Pipeline — Summary](#610-the-full-variant-calling-pipeline--summary)

---

## 6.1 What is a Variant?

A variant is any position in the genome where the sequenced sample differs from the reference genome. It is important to understand that the reference genome is not a "normal" or "correct" sequence — it is simply the agreed-upon standard. A variant simply means a difference, not necessarily a disease-causing change.

### Types of Variants

**Single Nucleotide Polymorphism (SNP)**  
A single base change at a specific position. The most common type of variant in the human genome — approximately 4–5 million SNPs distinguish any two unrelated people.

```
Reference:  ...ATCG[A]GCTA...
Sample:     ...ATCG[G]GCTA...
                    ↑ SNP: A→G
```

**Insertion / Deletion (InDel)**  
A small insertion or deletion of 1–50 bp. InDels in coding regions that are not multiples of 3 bp cause frameshifts, which can severely disrupt protein function.

```
Reference:  ...ATCG----GCTA...   (deletion)
Sample:     ...ATCGTTTTGCTA...
                 ↑↑↑↑ 4 bp insertion

Reference:  ...ATCGTTTTGCTA...   (deletion)
Sample:     ...ATCG----GCTA...
                 ↑↑↑↑ 4 bp deletion
```

**Copy Number Variant (CNV)**  
A segment of the genome (typically > 1 kb) that is duplicated or deleted. CNVs affect gene dosage — having 3 copies of a gene instead of the normal 2 can alter expression levels significantly.

**Structural Variant (SV)**  
Large genomic rearrangements (> 50 bp) including:
- Inversions — a segment is reversed in orientation
- Translocations — a segment moves to a different chromosome
- Large duplications or deletions

> **Workshop focus:** This module focuses on **SNPs and small InDels** using the GATK HaplotypeCaller — the most commonly performed variant calling task in research and clinical genomics. CNV and SV calling require different tools and are beyond the scope of this workshop.

### Variant Allele Frequency (VAF)

VAF is the proportion of reads at a position that carry the alternate allele:

```
VAF = (reads supporting ALT) / (total reads at position)
```

| VAF | Interpretation |
|---|---|
| ~0.5 (50%) | Heterozygous germline variant (one copy changed) |
| ~1.0 (100%) | Homozygous germline variant (both copies changed) |
| 0.01–0.49 | Somatic variant or mosaic variant |
| < 0.01 | Likely sequencing noise — requires very high coverage to call |

---

## 6.2 The VCF File Format

VCF (Variant Call Format) is the universal standard for storing variant calls. Understanding VCF is essential — every downstream analysis tool works with VCF files.

### VCF Structure

A VCF file has two sections:

**1. Header lines** — begin with `##`; describe the reference genome, tool versions, and define all INFO and FORMAT fields used in the file

**2. Data lines** — one line per variant position

### VCF Header

```
##fileformat=VCFv4.2
##GATK version=4.4.0.0
##reference=chr20.fa
##contig=<ID=chr20,length=64444167>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO                    FORMAT    SAMPLE1
```

### VCF Data Lines — The 8 Mandatory Fields

```
chr20  1234567  rs12345  A  G  892.77  PASS  DP=87;AF=0.51  GT:AD:DP:GQ  0/1:43,44:87:99
```

| Field | Example | Description |
|---|---|---|
| CHROM | chr20 | Chromosome |
| POS | 1234567 | 1-based position of the variant |
| ID | rs12345 | Known variant ID from dbSNP (`.` if unknown) |
| REF | A | Reference allele at this position |
| ALT | G | Alternate allele(s) observed in the sample |
| QUAL | 892.77 | Phred-scaled quality score for the variant call |
| FILTER | PASS | PASS if variant passed all filters; filter name if failed |
| INFO | DP=87;AF=0.51 | Variant-level annotations (depth, allele frequency, etc.) |

### The FORMAT and Sample Fields

After the 8 mandatory fields, VCF files contain a FORMAT field and one column per sample:

```
FORMAT:   GT:AD:DP:GQ
SAMPLE1:  0/1:43,44:87:99
```

| Tag | Value | Meaning |
|---|---|---|
| GT | 0/1 | Genotype — 0=reference allele, 1=first ALT allele |
| AD | 43,44 | Allelic depth — 43 reads support REF, 44 support ALT |
| DP | 87 | Total read depth at this position |
| GQ | 99 | Genotype quality — confidence in the GT call (Phred-scaled) |
| PL | 0,99,999 | Phred-scaled likelihoods for each possible genotype |

### Common Genotype Codes

| GT | Meaning |
|---|---|
| 0/0 | Homozygous reference — same as reference genome |
| 0/1 | Heterozygous — one reference, one alternate allele |
| 1/1 | Homozygous alternate — both alleles differ from reference |
| 0/2 | Heterozygous with second alternate allele |
| ./. | Missing genotype — insufficient coverage to call |

### Multi-allelic Sites

Some positions have more than one alternate allele:
```
chr20  5678901  .  A  G,T  .  PASS  .  GT  1/2
```
Here the sample is heterozygous for two different alternate alleles (G and T) — neither copy matches the reference.

---

## 6.3 Germline vs Somatic Variant Calling

The calling strategy differs fundamentally depending on whether you are looking for germline or somatic variants.

### Germline Variant Calling

Germline variants are present in every cell of the body — they were inherited from parents or arose in the germline before fertilisation. They are present at approximately 50% VAF (heterozygous) or 100% VAF (homozygous).

**Characteristics:**
- High VAF (0.5 or 1.0)
- Present in all tissues of the individual
- Relevant for inherited disease, pharmacogenomics, population genetics
- Requires one sample (tumour/normal pairing not needed)

**Tool:** GATK HaplotypeCaller (covered in detail in this module)

### Somatic Variant Calling

Somatic variants arise in individual cells during a person's lifetime — they are present in only a subset of cells (the tumour clone). In cancer genomics, the goal is to find mutations present in tumour cells but absent from normal (germline) cells.

**Characteristics:**
- Low VAF (0.05–0.40 typically) — only present in tumour cells
- Require a matched normal sample for comparison (tumour/normal pair)
- Relevant for cancer diagnosis, treatment selection, monitoring
- Harder to call — require higher coverage (60–100×) and sensitive callers

**Tool:** GATK Mutect2 (mentioned here but not the focus of this module)

| | Germline | Somatic |
|---|---|---|
| VAF | ~0.5 or ~1.0 | 0.01–0.49 |
| Coverage needed | 30× | 60–100× |
| Samples needed | 1 (proband) | 2 (tumour + normal) |
| Primary tool | HaplotypeCaller | Mutect2 |
| Application | Inherited disease, GWAS | Cancer genomics |

---

## 6.4 The GATK HaplotypeCaller — How it Works

### Overview

The GATK HaplotypeCaller (HC) is the gold-standard tool for germline SNP and indel calling. Rather than simply comparing each base in each read to the reference, it uses a more sophisticated approach based on local de novo assembly and haplotype-based genotyping.

### Step 1 — Active Region Detection

HaplotypeCaller scans the BAM file looking for **active regions** — windows of the genome where the reads show evidence of variation (mismatches, insertions, deletions, or soft clips above a threshold). Regions where reads match the reference perfectly are skipped entirely, saving computation.

### Step 2 — Local De Novo Assembly

Within each active region, HaplotypeCaller performs **local de novo assembly** using a De Bruijn graph. Rather than relying on the alignment produced by BWA-MEM2, it reconstructs all possible haplotypes (sequences) supported by the reads in that region. This is the key innovation of HaplotypeCaller — it can discover variants that short-read aligners miss or misrepresent, particularly complex indels and variants in regions where multiple nearby variants occur on the same haplotype.

### Step 3 — Pairwise Alignment of Haplotypes to Reference

The assembled haplotypes are aligned back to the reference using Smith-Waterman alignment. Each possible haplotype is evaluated — this is where SNPs and indels are identified relative to the reference.

### Step 4 — Likelihood Calculation and Genotyping

Each read is assigned a likelihood of having come from each haplotype. Using these likelihoods, HaplotypeCaller calculates the probability of each possible genotype (0/0, 0/1, 1/1) at each variant site using Bayesian statistics. The genotype with the highest posterior probability is assigned.

### Running HaplotypeCaller

**Single-sample mode (produces VCF directly):**

```bash
gatk HaplotypeCaller \
    -R reference/chr20.fa \
    -I aligned/sample_bqsr.bam \
    -O variants/sample_raw.vcf.gz \
    --dbsnp reference/dbsnp138_chr20.vcf.gz \
    2> logs/haplotypecaller.log
```

**GVCF mode (recommended for multi-sample projects):**

In GVCF mode, HaplotypeCaller produces a Genomic VCF that contains evidence for every position in the genome — not just variant sites. This allows multiple samples to be jointly genotyped later, which improves sensitivity and accuracy.

```bash
gatk HaplotypeCaller \
    -R reference/chr20.fa \
    -I aligned/sample_bqsr.bam \
    -O variants/sample.g.vcf.gz \
    -ERC GVCF \
    --dbsnp reference/dbsnp138_chr20.vcf.gz \
    2> logs/haplotypecaller_gvcf.log
```

**Joint genotyping from multiple GVCFs:**

```bash
# Combine GVCFs from multiple samples
gatk CombineGVCFs \
    -R reference/chr20.fa \
    -V variants/sample1.g.vcf.gz \
    -V variants/sample2.g.vcf.gz \
    -V variants/sample3.g.vcf.gz \
    -O variants/combined.g.vcf.gz

# Joint genotyping
gatk GenotypeGVCFs \
    -R reference/chr20.fa \
    -V variants/combined.g.vcf.gz \
    -O variants/joint_genotyped.vcf.gz \
    --dbsnp reference/dbsnp138_chr20.vcf.gz
```

---

## 6.5 Variant Filtering — VQSR and Hard Filters

Raw variant calls from HaplotypeCaller contain both true variants and false positives — sequencing errors, alignment artefacts, and PCR errors that the caller could not distinguish from real variants. Filtering is essential before any downstream analysis.

### Method 1 — Variant Quality Score Recalibration (VQSR)

VQSR is the recommended filtering method for large datasets (> 30 WGS samples or > 50 WES samples). It uses machine learning to build a model of what true variants look like based on a set of known truth variants (from HapMap, 1000 Genomes, Omni), then scores each variant by how well it matches the model.

VQSR requires a large number of variants to train the model — it is **not suitable for single samples or small panels**. For the workshop, we use hard filters instead.

```bash
# SNP VQSR — Step 1: Build recalibration model
gatk VariantRecalibrator \
    -R reference/chr20.fa \
    -V variants/joint_genotyped.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_snps.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp138.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O variants/snps_recal.vcf.gz \
    --tranches-file variants/snps.tranches

# SNP VQSR — Step 2: Apply recalibration
gatk ApplyVQSR \
    -R reference/chr20.fa \
    -V variants/joint_genotyped.vcf.gz \
    --recal-file variants/snps_recal.vcf.gz \
    --tranches-file variants/snps.tranches \
    --truth-sensitivity-filter-level 99.5 \
    -mode SNP \
    -O variants/vqsr_snps_filtered.vcf.gz
```

### Method 2 — Hard Filters (for small datasets and this workshop)

For single samples or small cohorts, GATK recommends applying hard filters — simple threshold-based filters on variant quality metrics. These are less powerful than VQSR but appropriate when there are insufficient variants to train a recalibration model.

**First, split SNPs and indels — they require different filter thresholds:**

```bash
# Extract SNPs only
gatk SelectVariants \
    -R reference/chr20.fa \
    -V variants/sample_raw.vcf.gz \
    --select-type-to-include SNP \
    -O variants/sample_snps_raw.vcf.gz

# Extract indels only
gatk SelectVariants \
    -R reference/chr20.fa \
    -V variants/sample_raw.vcf.gz \
    --select-type-to-include INDEL \
    -O variants/sample_indels_raw.vcf.gz
```

**Filter SNPs:**

```bash
gatk VariantFiltration \
    -R reference/chr20.fa \
    -V variants/sample_snps_raw.vcf.gz \
    --filter-expression "QD < 2.0"      --filter-name "QD2" \
    --filter-expression "FS > 60.0"     --filter-name "FS60" \
    --filter-expression "MQ < 40.0"     --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0"     --filter-name "SOR3" \
    -O variants/sample_snps_filtered.vcf.gz
```

**Filter indels:**

```bash
gatk VariantFiltration \
    -R reference/chr20.fa \
    -V variants/sample_indels_raw.vcf.gz \
    --filter-expression "QD < 2.0"      --filter-name "QD2" \
    --filter-expression "FS > 200.0"    --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0"    --filter-name "SOR10" \
    -O variants/sample_indels_filtered.vcf.gz
```

### Understanding the Filter Metrics

| Metric | Full Name | What it measures | SNP threshold | InDel threshold |
|---|---|---|---|---|
| QD | QualByDepth | Variant quality normalised by depth — low QD suggests a weak call | < 2.0 | < 2.0 |
| FS | FisherStrand | Strand bias — true variants appear on both strands equally | > 60.0 | > 200.0 |
| MQ | RMSMappingQuality | Average mapping quality of reads supporting the variant | < 40.0 | N/A |
| MQRankSum | MappingQualityRankSumTest | Difference in MAPQ between REF and ALT reads | < -12.5 | N/A |
| ReadPosRankSum | ReadPosRankSumTest | Whether ALT allele reads cluster at read ends (artefact sign) | < -8.0 | < -20.0 |
| SOR | StrandOddsRatio | More robust strand bias metric than FS | > 3.0 | > 10.0 |

### Merging Filtered SNPs and InDels

```bash
gatk MergeVcfs \
    -I variants/sample_snps_filtered.vcf.gz \
    -I variants/sample_indels_filtered.vcf.gz \
    -O variants/sample_filtered.vcf.gz
```

### Extracting PASS Variants Only

```bash
# Keep only variants that passed all filters
gatk SelectVariants \
    -R reference/chr20.fa \
    -V variants/sample_filtered.vcf.gz \
    --exclude-filtered \
    -O variants/sample_PASS.vcf.gz

# Count variants
echo "Total raw variants:"
bcftools stats variants/sample_raw.vcf.gz | grep "^SN" | grep "number of records"

echo "PASS variants:"
bcftools stats variants/sample_PASS.vcf.gz | grep "^SN" | grep "number of records"
```

---

## 6.6 Variant Annotation

A filtered VCF file tells you *where* variants are in the genome but not *what they mean* biologically. Variant annotation adds functional context — is this variant in a gene? Does it change a protein? Is it known in the population? Is it in ClinVar as pathogenic?

### SnpEff — Functional Effect Prediction

SnpEff predicts the functional consequence of each variant based on gene models. It assigns each variant an effect such as:

| Effect | Meaning |
|---|---|
| missense_variant | Changes one amino acid to another |
| synonymous_variant | Changes codon but not amino acid (silent) |
| stop_gained | Introduces a premature stop codon |
| frameshift_variant | InDel not divisible by 3 — disrupts reading frame |
| splice_region_variant | Near a splice site — may affect splicing |
| intron_variant | Within an intron — usually benign |
| upstream_gene_variant | Upstream of a gene — may affect regulation |

```bash
# Install SnpEff
conda install -c bioconda snpeff -y

# Download GRCh38 database
snpEff download GRCh38.99

# Annotate variants
snpEff ann \
    -v GRCh38.99 \
    -stats reports/sample_snpeff_summary.html \
    variants/sample_PASS.vcf.gz \
    > variants/sample_snpeff.vcf

# Compress and index
bgzip variants/sample_snpeff.vcf
tabix -p vcf variants/sample_snpeff.vcf.gz
```

### ANNOVAR — Population Frequency and Database Annotation

ANNOVAR annotates variants with information from population databases — giving you allele frequencies in different populations, ClinVar classification, COSMIC cancer entries, and functional prediction scores.

```bash
# Download ANNOVAR databases (requires registration at annovar website)
# humandb/ contains downloaded annotation databases

# Run ANNOVAR table annotation
table_annovar.pl \
    variants/sample_PASS.vcf.gz humandb/ \
    -buildver hg38 \
    -out variants/sample_annovar \
    -remove \
    -protocol refGene,cytoBand,gnomad30_genome,clinvar_20221231,cosmic97_coding \
    -operation g,r,f,f,f \
    -nastring . \
    -vcfinput \
    -polish
```

**Key databases used:**

| Database | What it adds |
|---|---|
| refGene | Gene name, transcript, amino acid change |
| gnomAD | Population allele frequencies across 125,000+ genomes |
| ClinVar | Clinical significance (pathogenic, benign, VUS) |
| COSMIC | Cancer somatic mutation database |
| CADD | Combined Annotation Dependent Depletion score — overall deleteriousness prediction |
| dbSNP | Known variant rsID |

### Ensembl VEP (Variant Effect Predictor)

VEP is the Ensembl alternative to ANNOVAR — widely used in clinical and research settings, with a web interface as well as a command-line tool.

```bash
# Install VEP
conda install -c bioconda ensembl-vep -y

# Download cache
vep_install -a cf -s homo_sapiens -y GRCh38

# Run VEP
vep \
    --input_file variants/sample_PASS.vcf.gz \
    --output_file variants/sample_vep.vcf \
    --vcf \
    --everything \
    --fork 4 \
    --cache \
    --offline \
    --assembly GRCh38 \
    --stats_file reports/sample_vep_summary.html
```

`--everything` enables all VEP annotation plugins including SIFT, PolyPhen-2, regulatory annotations, and population frequencies.

---

## 6.7 Interpreting Annotated Variants

### Variant Prioritisation

A typical 30× WGS analysis of a human sample produces approximately **4–5 million variants** relative to the reference genome. The vast majority are common polymorphisms — harmless differences present in many people. The challenge in clinical genomics is narrowing this down to the small number of variants that may be clinically or biologically relevant.

A standard prioritisation funnel for rare disease genomics:

```
~4,500,000 total variants
        ↓  Remove common variants (gnomAD AF > 1%)
~100,000 rare variants
        ↓  Keep only coding + splice site variants
~10,000 coding rare variants
        ↓  Keep only high-impact (stop, frameshift, missense)
~2,000 high-impact rare variants
        ↓  Filter by inheritance pattern (de novo, recessive, dominant)
~50–200 candidate variants
        ↓  Manual review, literature search, functional evidence
~1–10 likely pathogenic variants
```

### Key Annotation Fields to Examine

**gnomAD allele frequency (AF)**  
If a variant is present at > 1% frequency in gnomAD, it is unlikely to cause a rare Mendelian disease. Common variants (AF > 5%) are almost certainly benign in the context of rare disease.

```
gnomAD_genome_ALL = 0.0001   → rare variant, investigate further
gnomAD_genome_ALL = 0.35     → common polymorphism, likely benign
gnomAD_genome_ALL = .        → not in gnomAD — potentially novel
```

**ClinVar classification**  
ClinVar aggregates clinical variant interpretations from laboratories worldwide:

| ClinVar classification | Meaning |
|---|---|
| Pathogenic | Strong evidence this variant causes disease |
| Likely pathogenic | Good evidence this variant causes disease |
| Variant of Uncertain Significance (VUS) | Insufficient evidence to classify |
| Likely benign | Good evidence this variant does not cause disease |
| Benign | Strong evidence this variant does not cause disease |

**CADD score**  
CADD (Combined Annotation Dependent Depletion) integrates many annotations into a single score predicting variant deleteriousness. Scores are Phred-scaled:

| CADD score | Interpretation |
|---|---|
| < 10 | Bottom 90% of variants — likely benign |
| ≥ 10 | Top 10% — potentially functional |
| ≥ 20 | Top 1% — likely deleterious |
| ≥ 30 | Top 0.1% — highly deleterious |

**SIFT and PolyPhen-2**  
These tools predict whether a missense variant (amino acid change) disrupts protein function:
- SIFT: `deleterious` (score < 0.05) or `tolerated` (score ≥ 0.05)
- PolyPhen-2: `probably_damaging`, `possibly_damaging`, or `benign`

> ⚠️ **Important:** No single annotation tool is definitive. A variant classified as "probably damaging" by PolyPhen-2 and "deleterious" by SIFT is a candidate for further investigation — not a confirmed pathogenic variant. Always integrate multiple lines of evidence and consult the primary literature before drawing clinical conclusions.

### Parsing Annotated VCF Files with bcftools

```bash
# Install bcftools
conda install -c bioconda bcftools -y

# View all PASS variants
bcftools view -f PASS variants/sample_filtered.vcf.gz | less

# Count SNPs and indels separately
bcftools stats variants/sample_PASS.vcf.gz | grep "^SN"

# Extract only coding variants (missense, stop-gain, frameshift)
bcftools view variants/sample_snpeff.vcf.gz | \
    grep -E "missense_variant|stop_gained|frameshift_variant" | \
    head -20

# Filter by allele frequency — keep only rare variants (AF < 0.01)
bcftools view -i 'INFO/AF < 0.01' variants/sample_PASS.vcf.gz \
    -o variants/sample_rare.vcf.gz -Oz

# Convert VCF to tab-delimited table for easy viewing in spreadsheet
gatk VariantsToTable \
    -V variants/sample_snpeff.vcf.gz \
    -F CHROM -F POS -F REF -F ALT -F FILTER -F QUAL \
    -GF GT -GF AD -GF DP -GF GQ \
    -O variants/sample_variants_table.tsv
```

---

## 6.8 Hands-on Exercises

### Setup

```bash
mkdir -p ~/ngs_workshop/module6/{variants,reports,logs}
cd ~/ngs_workshop/module6

# Install tools
conda install -c bioconda gatk4 snpeff bcftools -y

# Symlink the BAM from Module 5
ln -s ~/ngs_workshop/module5/aligned ./aligned
ln -s ~/ngs_workshop/module5/reference ./reference
```

---

### Exercise 6.1 — Run HaplotypeCaller

```bash
cd ~/ngs_workshop/module6

gatk HaplotypeCaller \
    -R reference/chr20.fa \
    -I aligned/sample_bqsr.bam \
    -O variants/sample_raw.vcf.gz \
    --dbsnp reference/dbsnp138_chr20.vcf.gz \
    -L chr20 \
    2> logs/haplotypecaller.log

echo "Variant calling complete"
echo "Raw variants called:"
bcftools stats variants/sample_raw.vcf.gz | grep "number of records"
```

**Questions:**
1. How many raw variant sites were called?
2. How long did HaplotypeCaller take?
3. Look at the log — how many active regions were processed?

---

### Exercise 6.2 — Inspect the Raw VCF

```bash
# View the VCF header
bcftools view -h variants/sample_raw.vcf.gz | tail -20

# View the first 10 variant records
bcftools view variants/sample_raw.vcf.gz | grep -v "^#" | head -10

# Count SNPs vs indels
echo "SNPs:"
bcftools stats variants/sample_raw.vcf.gz | grep "number of SNPs"
echo "Indels:"
bcftools stats variants/sample_raw.vcf.gz | grep "number of indels"

# Check the Ti/Tv ratio (transitions to transversions)
# For WGS: expected ~2.0-2.1; for WES: ~2.5-3.0
# Values far outside this range suggest poor quality calls
bcftools stats variants/sample_raw.vcf.gz | grep "Ts/Tv"
```

**What is Ti/Tv ratio?**  
Transitions (Ti) are purine↔purine or pyrimidine↔pyrimidine changes (A↔G, C↔T). Transversions (Tv) are purine↔pyrimidine changes. In human genomes, transitions are biologically more common, giving a Ti/Tv of ~2.1 for WGS. A ratio well below 2.0 suggests an excess of false positive calls (random errors have a Ti/Tv of ~0.5).

---

### Exercise 6.3 — Apply Hard Filters

```bash
# Split SNPs and indels
gatk SelectVariants -R reference/chr20.fa \
    -V variants/sample_raw.vcf.gz \
    --select-type-to-include SNP \
    -O variants/sample_snps_raw.vcf.gz

gatk SelectVariants -R reference/chr20.fa \
    -V variants/sample_raw.vcf.gz \
    --select-type-to-include INDEL \
    -O variants/sample_indels_raw.vcf.gz

# Filter SNPs
gatk VariantFiltration -R reference/chr20.fa \
    -V variants/sample_snps_raw.vcf.gz \
    --filter-expression "QD < 2.0"          --filter-name "QD2" \
    --filter-expression "FS > 60.0"          --filter-name "FS60" \
    --filter-expression "MQ < 40.0"          --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"  --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0"          --filter-name "SOR3" \
    -O variants/sample_snps_filtered.vcf.gz

# Filter indels
gatk VariantFiltration -R reference/chr20.fa \
    -V variants/sample_indels_raw.vcf.gz \
    --filter-expression "QD < 2.0"           --filter-name "QD2" \
    --filter-expression "FS > 200.0"         --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0"         --filter-name "SOR10" \
    -O variants/sample_indels_filtered.vcf.gz

# Merge back together
gatk MergeVcfs \
    -I variants/sample_snps_filtered.vcf.gz \
    -I variants/sample_indels_filtered.vcf.gz \
    -O variants/sample_filtered.vcf.gz

# Extract PASS only
gatk SelectVariants -R reference/chr20.fa \
    -V variants/sample_filtered.vcf.gz \
    --exclude-filtered \
    -O variants/sample_PASS.vcf.gz

echo "Variants before filtering:"
bcftools stats variants/sample_raw.vcf.gz | grep "number of records"
echo "Variants after filtering (PASS only):"
bcftools stats variants/sample_PASS.vcf.gz | grep "number of records"
```

**Questions:**
1. How many variants were removed by filtering?
2. What percentage of raw variants passed all filters?
3. Did the Ti/Tv ratio improve after filtering?

---

### Exercise 6.4 — Annotate with SnpEff

```bash
# Download GRCh38 database if not already present
snpEff download GRCh38.99

# Annotate
snpEff ann \
    -v GRCh38.99 \
    -stats reports/snpeff_summary.html \
    variants/sample_PASS.vcf.gz \
    > variants/sample_snpeff.vcf

bgzip variants/sample_snpeff.vcf
tabix -p vcf variants/sample_snpeff.vcf.gz

# Count variants by effect
echo "Missense variants:"
bcftools view variants/sample_snpeff.vcf.gz | grep "missense_variant" | wc -l

echo "Stop-gained variants:"
bcftools view variants/sample_snpeff.vcf.gz | grep "stop_gained" | wc -l

echo "Frameshift variants:"
bcftools view variants/sample_snpeff.vcf.gz | grep "frameshift_variant" | wc -l

# Open the SnpEff HTML summary
xdg-open reports/snpeff_summary.html
```

---

### Exercise 6.5 — Convert to Table and Explore

```bash
# Convert annotated VCF to a TSV table
gatk VariantsToTable \
    -V variants/sample_snpeff.vcf.gz \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER \
    -GF GT -GF AD -GF DP -GF GQ \
    -O variants/sample_variants_table.tsv

# View the table
head -5 variants/sample_variants_table.tsv | column -t

# Count heterozygous vs homozygous variants
echo "Heterozygous variants (0/1):"
grep -c "0/1" variants/sample_variants_table.tsv

echo "Homozygous alt variants (1/1):"
grep -c "1/1" variants/sample_variants_table.tsv
```

---

### Exercise 6.6 — Write the Full Variant Calling Script

```bash
cat > ~/ngs_workshop/module6/run_variant_calling.sh << 'EOF'
#!/bin/bash
# Full variant calling pipeline: BAM → annotated VCF
# Usage: bash run_variant_calling.sh <sample_name>

set -euo pipefail

SAMPLE=$1
REF=reference/chr20.fa
DBSNP=reference/dbsnp138_chr20.vcf.gz
THREADS=4

mkdir -p variants reports logs

echo "[$(date)] Starting variant calling for: ${SAMPLE}"

# Step 1 — HaplotypeCaller
echo "[$(date)] Step 1/5: Running HaplotypeCaller..."
gatk HaplotypeCaller \
    -R ${REF} \
    -I aligned/${SAMPLE}_bqsr.bam \
    -O variants/${SAMPLE}_raw.vcf.gz \
    --dbsnp ${DBSNP} \
    -L chr20 \
    2> logs/${SAMPLE}_haplotypecaller.log

# Step 2 — Split, filter, merge
echo "[$(date)] Step 2/5: Filtering variants..."
gatk SelectVariants -R ${REF} -V variants/${SAMPLE}_raw.vcf.gz \
    --select-type-to-include SNP -O variants/${SAMPLE}_snps_raw.vcf.gz
gatk SelectVariants -R ${REF} -V variants/${SAMPLE}_raw.vcf.gz \
    --select-type-to-include INDEL -O variants/${SAMPLE}_indels_raw.vcf.gz

gatk VariantFiltration -R ${REF} -V variants/${SAMPLE}_snps_raw.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    -O variants/${SAMPLE}_snps_filtered.vcf.gz

gatk VariantFiltration -R ${REF} -V variants/${SAMPLE}_indels_raw.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10" \
    -O variants/${SAMPLE}_indels_filtered.vcf.gz

gatk MergeVcfs \
    -I variants/${SAMPLE}_snps_filtered.vcf.gz \
    -I variants/${SAMPLE}_indels_filtered.vcf.gz \
    -O variants/${SAMPLE}_filtered.vcf.gz

gatk SelectVariants -R ${REF} -V variants/${SAMPLE}_filtered.vcf.gz \
    --exclude-filtered -O variants/${SAMPLE}_PASS.vcf.gz

# Step 3 — Annotate with SnpEff
echo "[$(date)] Step 3/5: Annotating with SnpEff..."
snpEff ann -v GRCh38.99 \
    -stats reports/${SAMPLE}_snpeff_summary.html \
    variants/${SAMPLE}_PASS.vcf.gz \
    > variants/${SAMPLE}_snpeff.vcf
bgzip variants/${SAMPLE}_snpeff.vcf
tabix -p vcf variants/${SAMPLE}_snpeff.vcf.gz

# Step 4 — Convert to table
echo "[$(date)] Step 4/5: Converting to table..."
gatk VariantsToTable \
    -V variants/${SAMPLE}_snpeff.vcf.gz \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER \
    -GF GT -GF AD -GF DP -GF GQ \
    -O variants/${SAMPLE}_variants_table.tsv

# Step 5 — Summary stats
echo "[$(date)] Step 5/5: Summary statistics..."
echo "" > reports/${SAMPLE}_variant_summary.txt
echo "=== Variant Summary: ${SAMPLE} ===" >> reports/${SAMPLE}_variant_summary.txt
echo "Raw variants:" >> reports/${SAMPLE}_variant_summary.txt
bcftools stats variants/${SAMPLE}_raw.vcf.gz | grep "number of records" >> reports/${SAMPLE}_variant_summary.txt
echo "PASS variants:" >> reports/${SAMPLE}_variant_summary.txt
bcftools stats variants/${SAMPLE}_PASS.vcf.gz | grep "number of records" >> reports/${SAMPLE}_variant_summary.txt
echo "Ti/Tv ratio (PASS):" >> reports/${SAMPLE}_variant_summary.txt
bcftools stats variants/${SAMPLE}_PASS.vcf.gz | grep "Ts/Tv" >> reports/${SAMPLE}_variant_summary.txt

cat reports/${SAMPLE}_variant_summary.txt

echo "[$(date)] Pipeline complete!"
echo "Annotated VCF: variants/${SAMPLE}_snpeff.vcf.gz"
echo "Variant table: variants/${SAMPLE}_variants_table.tsv"
EOF

chmod +x ~/ngs_workshop/module6/run_variant_calling.sh

# Run it
bash ~/ngs_workshop/module6/run_variant_calling.sh SRR6821753
```

---

## 6.9 Common Problems and Troubleshooting

| Problem | Likely Cause | Solution |
|---|---|---|
| GATK error: `no read groups found` | BAM missing `@RG` header | Re-run alignment with `-R` read group tag (Module 5) |
| Very few variants called (< 1000 for chr20) | Wrong reference genome; low coverage; only chr20 data | Confirm reference matches sample; check coverage |
| Ti/Tv ratio < 1.5 | Too many false positives | Tighten hard filter thresholds; check alignment quality |
| Ti/Tv ratio > 3.5 | Over-filtering or WES data | Expected for WES; for WGS check filter thresholds |
| All variants filtered out | Filter thresholds too strict | Loosen one filter at a time and check impact |
| SnpEff `database not found` | Database not downloaded | Run `snpEff download GRCh38.99` |
| VCF indexing error | VCF not bgzipped | Use `bgzip` not `gzip`; then `tabix -p vcf file.vcf.gz` |
| GATK OutOfMemory error | Java heap too small | Add `-Xmx8g` flag: `gatk --java-options "-Xmx8g" HaplotypeCaller` |

---

## 6.10 The Full Variant Calling Pipeline — Summary

```
Analysis-ready BAM (from Module 5)
        ↓
1. GATK HaplotypeCaller
   → raw VCF (all candidate variants)
        ↓
2. SelectVariants (split SNPs / indels)
        ↓
3. VariantFiltration (hard filters)
        ↓
4. MergeVcfs + SelectVariants PASS
   → filtered VCF (high-confidence variants)
        ↓
5. SnpEff / VEP / ANNOVAR annotation
   → annotated VCF (variants with biological context)
        ↓
6. Prioritisation and interpretation
   → candidate variants for further investigation
        ↓
      Module 7: Visualisation and Interpretation
```

---

> 💡 **Key takeaway for this module:** Variant calling is not a single step — it is a pipeline of calling, filtering, and annotation, each of which requires careful parameter choices. The Ti/Tv ratio is your best quick sanity check after filtering. An annotated VCF file is a starting point for biological interpretation, not an endpoint — always validate important findings with orthogonal methods (Sanger sequencing, functional assays) before drawing conclusions.

---

**Previous:** [Module 5 — Sequence Alignment and Mapping](./module5_alignment.md)  
**Next:** [Module 7 — Visualization and Interpretation](./module7_visualization.md)
