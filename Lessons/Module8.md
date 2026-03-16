# Module 8: Advanced NGS Pipelines and Workflow Management

> **NGS Workshop — 2-Day Intensive Program for UG & PG Students**  
> Day 2 | Module 8 of 4 | Final Module

> **Context:** You have now run every step of an NGS pipeline manually — QC, preprocessing, alignment, variant calling, annotation, and visualisation. This module teaches you how to automate, scale, and share those pipelines professionally using workflow management tools.

---

## Table of Contents

- [8.1 The Problem with Manual Pipelines](#81-the-problem-with-manual-pipelines)
- [8.2 What is a Workflow Manager?](#82-what-is-a-workflow-manager)
- [8.3 Snakemake](#83-snakemake)
- [8.4 Nextflow](#84-nextflow)
- [8.5 nf-core — Community NGS Pipelines](#85-nf-core--community-ngs-pipelines)
- [8.6 Containerisation — Docker and Singularity](#86-containerisation--docker-and-singularity)
- [8.7 Running Pipelines on HPC and Cloud](#87-running-pipelines-on-hpc-and-cloud)
- [8.8 Pipeline Best Practices](#88-pipeline-best-practices)
- [8.9 Hands-on Exercises](#89-hands-on-exercises)
- [8.10 Where to Go From Here](#810-where-to-go-from-here)

---

## 8.1 The Problem with Manual Pipelines

Over the past two days, you have run each step of an NGS pipeline manually — one command at a time, for one sample. In research and clinical practice, this approach breaks down quickly:

**Scale** — A typical research project involves 10–200 samples. Running each step manually for every sample is not feasible. A single RNA-seq study comparing two conditions might have 6 samples per group — 12 samples × 8 pipeline steps = 96 individual commands to run, monitor, and track.

**Reproducibility** — If you run commands manually, it is easy to make subtle mistakes — using slightly different parameters for different samples, accidentally overwriting files, or forgetting which version of a tool you used. When reviewers or collaborators ask you to reproduce your results six months later, can you?

**Error handling** — If one step fails (e.g. a disk is full, a node crashes on an HPC), a manual pipeline just stops. A workflow manager detects the failure, reports exactly which step and sample failed, and allows you to resume from where it left off without rerunning completed steps.

**Portability** — A pipeline that runs on your laptop should also run on an HPC cluster, a cloud instance, or a colleague's machine without rewriting all the commands.

Workflow managers solve all of these problems. They allow you to define your pipeline once — as a series of rules or processes — and then execute it across any number of samples, on any computing infrastructure, with full logging, error handling, and reproducibility.

---

## 8.2 What is a Workflow Manager?

A workflow manager is a tool that:

1. **Defines dependencies** between pipeline steps — it knows that alignment must happen before duplicate marking, which must happen before variant calling
2. **Executes steps in parallel** where possible — if you have 20 samples, it can align all 20 simultaneously rather than one at a time
3. **Tracks completed steps** — if a run fails halfway through, it resumes from the point of failure rather than starting over
4. **Manages compute resources** — allocates CPUs, memory, and time appropriately for each step
5. **Produces a full audit trail** — every command run, every tool version, every parameter is logged

### The Two Dominant Tools in Bioinformatics

| | Snakemake | Nextflow |
|---|---|---|
| **Language** | Python-based rules | Groovy/DSL2 |
| **Learning curve** | Easier — good for Python users | Steeper — but very powerful |
| **Best for** | Academic labs, smaller projects | Large-scale production pipelines |
| **Community pipelines** | Snakemake Workflow Catalog | nf-core (gold standard) |
| **Cloud/HPC support** | Yes | Yes — excellent |
| **Containerisation** | Conda, Singularity, Docker | Docker, Singularity, Conda |
| **India usage** | Very common in academia | Growing, especially in industry |

Both are excellent tools. This module covers both, with enough depth to get you started and enough context to choose the right one for your project.

---

## 8.3 Snakemake

### What is Snakemake?

Snakemake was developed by Johannes Köster and is written in Python. A Snakemake pipeline is defined in a `Snakefile` — a text file containing a series of **rules**, each describing one step of the pipeline (its input files, output files, and the shell command to run).

The key concept is **target-based execution** — you tell Snakemake what output you want, and it works backwards through the rules to figure out what needs to be run, in what order, and for which samples.

### Installation

```bash
conda install -c bioconda -c conda-forge snakemake -y

# Verify
snakemake --version
```

### Core Concepts

**Rules** — the building blocks of a Snakemake pipeline. Each rule has:
- `input` — the file(s) this step needs
- `output` — the file(s) this step produces
- `shell` or `run` — the command to execute

**Wildcards** — placeholders like `{sample}` that Snakemake replaces with actual sample names, allowing one rule to handle many samples.

**The `all` rule** — a special rule at the top of the Snakefile that lists all final output files. Snakemake uses this to determine what needs to be run.

### A Complete NGS Snakemake Pipeline

This Snakefile implements the full variant calling pipeline from Modules 3–6:

```python
# Snakefile — NGS Variant Calling Pipeline
# Run with: snakemake --cores 8 --use-conda

# ─────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────
SAMPLES = ["SRR6821753", "SRR6821754"]  # add all sample names here
REF     = "reference/chr20.fa"
DBSNP   = "reference/dbsnp138_chr20.vcf.gz"

# ─────────────────────────────────────────
# Rule all — defines the final targets
# ─────────────────────────────────────────
rule all:
    input:
        expand("results/variants/{sample}_snpeff.vcf.gz", sample=SAMPLES),
        expand("results/reports/{sample}_flagstat.txt",   sample=SAMPLES),
        "results/multiqc/multiqc_report.html"

# ─────────────────────────────────────────
# Step 1 — FastQC on raw reads
# ─────────────────────────────────────────
rule fastqc_raw:
    input:
        r1 = "raw_data/{sample}_R1.fastq.gz",
        r2 = "raw_data/{sample}_R2.fastq.gz"
    output:
        html_r1 = "results/qc/fastqc/{sample}_R1_fastqc.html",
        html_r2 = "results/qc/fastqc/{sample}_R2_fastqc.html",
        zip_r1  = "results/qc/fastqc/{sample}_R1_fastqc.zip",
        zip_r2  = "results/qc/fastqc/{sample}_R2_fastqc.zip"
    threads: 2
    shell:
        "fastqc {input.r1} {input.r2} -o results/qc/fastqc/ -t {threads}"

# ─────────────────────────────────────────
# Step 2 — Trimming with fastp
# ─────────────────────────────────────────
rule fastp:
    input:
        r1 = "raw_data/{sample}_R1.fastq.gz",
        r2 = "raw_data/{sample}_R2.fastq.gz"
    output:
        r1    = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2    = "results/trimmed/{sample}_R2_trimmed.fastq.gz",
        html  = "results/qc/fastp/{sample}_fastp.html",
        json  = "results/qc/fastp/{sample}_fastp.json"
    threads: 4
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --html {output.html} --json {output.json} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 36 \
            --trim_poly_g \
            --thread {threads}
        """

# ─────────────────────────────────────────
# Step 3 — Alignment with BWA-MEM2
# ─────────────────────────────────────────
rule align:
    input:
        r1  = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2  = "results/trimmed/{sample}_R2_trimmed.fastq.gz",
        ref = REF
    output:
        bam = "results/aligned/{sample}_sorted.bam",
        bai = "results/aligned/{sample}_sorted.bam.bai"
    threads: 8
    shell:
        """
        bwa-mem2 mem -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:lib1\\tPU:{wildcards.sample}.1" \
            {input.ref} {input.r1} {input.r2} \
            | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

# ─────────────────────────────────────────
# Step 4 — Mark duplicates
# ─────────────────────────────────────────
rule mark_duplicates:
    input:
        bam = "results/aligned/{sample}_sorted.bam"
    output:
        bam     = "results/aligned/{sample}_markdup.bam",
        bai     = "results/aligned/{sample}_markdup.bam.bai",
        metrics = "results/reports/{sample}_markdup_metrics.txt"
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
        samtools index {output.bam}
        """

# ─────────────────────────────────────────
# Step 5 — BQSR
# ─────────────────────────────────────────
rule bqsr:
    input:
        bam    = "results/aligned/{sample}_markdup.bam",
        ref    = REF,
        dbsnp  = DBSNP
    output:
        bam    = "results/aligned/{sample}_bqsr.bam",
        bai    = "results/aligned/{sample}_bqsr.bam.bai",
        table  = "results/reports/{sample}_recal.table"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} -R {input.ref} \
            --known-sites {input.dbsnp} \
            -O {output.table}
        gatk ApplyBQSR \
            -I {input.bam} -R {input.ref} \
            --bqsr-recal-file {output.table} \
            -O {output.bam}
        samtools index {output.bam}
        """

# ─────────────────────────────────────────
# Step 6 — Variant calling
# ─────────────────────────────────────────
rule haplotype_caller:
    input:
        bam   = "results/aligned/{sample}_bqsr.bam",
        ref   = REF,
        dbsnp = DBSNP
    output:
        vcf = "results/variants/{sample}_raw.vcf.gz"
    threads: 4
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --dbsnp {input.dbsnp} \
            -L chr20
        """

# ─────────────────────────────────────────
# Step 7 — Hard filter variants
# ─────────────────────────────────────────
rule filter_variants:
    input:
        vcf = "results/variants/{sample}_raw.vcf.gz",
        ref = REF
    output:
        vcf = "results/variants/{sample}_PASS.vcf.gz"
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} \
            --select-type-to-include SNP \
            -O results/variants/{wildcards.sample}_snps.vcf.gz

        gatk VariantFiltration -R {input.ref} \
            -V results/variants/{wildcards.sample}_snps.vcf.gz \
            --filter-expression "QD < 2.0"    --filter-name "QD2" \
            --filter-expression "FS > 60.0"   --filter-name "FS60" \
            --filter-expression "MQ < 40.0"   --filter-name "MQ40" \
            --filter-expression "SOR > 3.0"   --filter-name "SOR3" \
            -O results/variants/{wildcards.sample}_snps_filtered.vcf.gz

        gatk SelectVariants -R {input.ref} -V {input.vcf} \
            --select-type-to-include INDEL \
            -O results/variants/{wildcards.sample}_indels.vcf.gz

        gatk VariantFiltration -R {input.ref} \
            -V results/variants/{wildcards.sample}_indels.vcf.gz \
            --filter-expression "QD < 2.0"    --filter-name "QD2" \
            --filter-expression "FS > 200.0"  --filter-name "FS200" \
            --filter-expression "SOR > 10.0"  --filter-name "SOR10" \
            -O results/variants/{wildcards.sample}_indels_filtered.vcf.gz

        gatk MergeVcfs \
            -I results/variants/{wildcards.sample}_snps_filtered.vcf.gz \
            -I results/variants/{wildcards.sample}_indels_filtered.vcf.gz \
            -O results/variants/{wildcards.sample}_filtered.vcf.gz

        gatk SelectVariants -R {input.ref} \
            -V results/variants/{wildcards.sample}_filtered.vcf.gz \
            --exclude-filtered \
            -O {output.vcf}
        """

# ─────────────────────────────────────────
# Step 8 — Annotate with SnpEff
# ─────────────────────────────────────────
rule snpeff:
    input:
        vcf = "results/variants/{sample}_PASS.vcf.gz"
    output:
        vcf  = "results/variants/{sample}_snpeff.vcf.gz",
        html = "results/reports/{sample}_snpeff_summary.html"
    shell:
        """
        snpEff ann -v GRCh38.99 \
            -stats {output.html} \
            {input.vcf} \
            | bgzip > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# ─────────────────────────────────────────
# Step 9 — Alignment QC metrics
# ─────────────────────────────────────────
rule flagstat:
    input:
        bam = "results/aligned/{sample}_bqsr.bam"
    output:
        txt = "results/reports/{sample}_flagstat.txt"
    shell:
        "samtools flagstat {input.bam} > {output.txt}"

# ─────────────────────────────────────────
# Step 10 — MultiQC report
# ─────────────────────────────────────────
rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("results/qc/fastp/{sample}_fastp.json",     sample=SAMPLES),
        expand("results/reports/{sample}_flagstat.txt",    sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        "multiqc results/qc/ results/reports/ -o results/multiqc/ -f"
```

### Running the Pipeline

```bash
# Dry run first — shows what will be executed without running anything
snakemake --dry-run --cores 8

# Run with 8 cores
snakemake --cores 8

# Run with conda environments (each rule gets its own environment)
snakemake --cores 8 --use-conda

# Run only specific samples
snakemake --cores 8 results/variants/SRR6821753_snpeff.vcf.gz

# Generate a visualisation of the DAG (directed acyclic graph)
snakemake --dag | dot -Tpng > pipeline_dag.png
```

### Snakemake Rule with Conda Environment

Each rule can specify its own conda environment, ensuring tool version reproducibility:

```python
rule fastp:
    input:  ...
    output: ...
    conda:  "envs/fastp.yaml"
    shell:  "fastp ..."
```

```yaml
# envs/fastp.yaml
channels:
  - bioconda
  - conda-forge
dependencies:
  - fastp=0.23.4
```

---

## 8.4 Nextflow

### What is Nextflow?

Nextflow was developed by Paolo Di Tommaso at the Centre for Genomic Regulation (CRG), Barcelona. It uses a dataflow programming model — data flows between **processes** (pipeline steps) through **channels** (queues). Nextflow is designed from the ground up for scalability, and the same pipeline code runs without modification on a laptop, an HPC cluster, or AWS/Google Cloud.

Nextflow uses **DSL2** (Domain Specific Language 2), which allows pipelines to be modular — individual processes are defined separately and assembled into workflows, making them reusable across projects.

### Installation

```bash
# Using conda
conda install -c bioconda nextflow -y

# Or directly
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify
nextflow -version
```

### Core Concepts

**Process** — equivalent to a Snakemake rule; defines one step with inputs, outputs, and a script block.

**Channel** — a queue that passes data between processes; can be a list of files, sample pairs, or any data.

**Workflow** — assembles processes and channels into a complete pipeline.

**Executor** — defines where processes run: `local`, `slurm`, `awsbatch`, `google-lifesciences`, etc.

### A Nextflow Pipeline — Core Structure

```groovy
// main.nf — NGS Variant Calling Pipeline in Nextflow DSL2

nextflow.enable.dsl = 2

// ─────────────────────────────────────────
// Parameters
// ─────────────────────────────────────────
params.reads    = "raw_data/*_{R1,R2}.fastq.gz"
params.ref      = "reference/chr20.fa"
params.dbsnp    = "reference/dbsnp138_chr20.vcf.gz"
params.outdir   = "results"

// ─────────────────────────────────────────
// Process definitions
// ─────────────────────────────────────────

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${reads} -t 2
    """
}

process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.{html,json}"

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample_id}_R1_trimmed.fastq.gz \\
        -O ${sample_id}_R2_trimmed.fastq.gz \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        --detect_adapter_for_pe \\
        --length_required 36 \\
        --thread 4
    """
}

process ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam

    script:
    """
    bwa-mem2 mem -t 8 \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:lib1\\tPU:${sample_id}.1" \\
        ${ref} ${reads[0]} ${reads[1]} \\
        | samtools sort -@ 8 -o ${sample_id}_sorted.bam -
    samtools index ${sample_id}_sorted.bam
    """
}

process HAPLOTYPECALLER {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path dbsnp

    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz"), emit: vcf

    script:
    """
    gatk HaplotypeCaller \\
        -R ${ref} \\
        -I ${bam} \\
        -O ${sample_id}_raw.vcf.gz \\
        --dbsnp ${dbsnp} \\
        -L chr20
    """
}

// ─────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────
workflow {
    // Create channel of read pairs
    reads_ch = Channel
        .fromFilePairs(params.reads)
        .ifEmpty { error "No reads found matching: ${params.reads}" }

    ref_ch   = Channel.fromPath(params.ref)
    dbsnp_ch = Channel.fromPath(params.dbsnp)

    // Execute processes
    FASTQC(reads_ch)
    FASTP(reads_ch)
    ALIGN(FASTP.out.trimmed_reads, ref_ch)
    HAPLOTYPECALLER(ALIGN.out.bam, ref_ch, dbsnp_ch)
}
```

### Running a Nextflow Pipeline

```bash
# Run locally
nextflow run main.nf

# Run with a specific profile
nextflow run main.nf -profile conda

# Resume a failed run from the last successful step
nextflow run main.nf -resume

# Run with custom parameters
nextflow run main.nf \
    --reads "raw_data/*_{R1,R2}.fastq.gz" \
    --ref reference/chr20.fa \
    --outdir my_results

# View pipeline execution report
nextflow run main.nf -with-report report.html -with-timeline timeline.html
```

### Nextflow Configuration File

The `nextflow.config` file controls compute resources and execution profiles:

```groovy
// nextflow.config

profiles {

    // Run locally with conda
    conda {
        conda.enabled = true
        process.conda = "envs/ngs_tools.yaml"
    }

    // Run on SLURM HPC cluster
    slurm {
        process.executor = 'slurm'
        process.queue    = 'compute'
        process.memory   = '16 GB'
        process.cpus     = 8
        process.time     = '4h'
    }

    // Run on AWS Batch
    awsbatch {
        process.executor    = 'awsbatch'
        process.queue       = 'nextflow-queue'
        aws.region          = 'ap-south-1'   // Mumbai — good for India
        aws.batch.cliPath   = '/usr/local/bin/aws'
    }
}

// Resource defaults
process {
    withName: ALIGN           { cpus = 8;  memory = '16 GB'; time = '4h' }
    withName: HAPLOTYPECALLER { cpus = 4;  memory = '8 GB';  time = '6h' }
    withName: FASTP           { cpus = 4;  memory = '4 GB';  time = '1h' }
}
```

---

## 8.5 nf-core — Community NGS Pipelines

### What is nf-core?

nf-core (https://nf-co.re) is a community-curated collection of high-quality, peer-reviewed Nextflow pipelines for common NGS analyses. Rather than writing your own pipeline from scratch, you can use a battle-tested nf-core pipeline that implements best practices, handles edge cases, and has been validated on thousands of samples across dozens of institutions worldwide.

nf-core pipelines are:
- Fully documented with usage guides and parameter references
- Containerised — every tool runs in a Docker/Singularity container for perfect reproducibility
- Tested on every commit with continuous integration
- Regularly updated as new tools and best practices emerge

### Most Relevant nf-core Pipelines for This Workshop

| Pipeline | What it does | Relevant modules |
|---|---|---|
| **nf-core/sarek** | Germline and somatic variant calling (WGS/WES) | Modules 5, 6 |
| **nf-core/rnaseq** | RNA-seq QC, alignment, and quantification | — |
| **nf-core/chipseq** | ChIP-seq peak calling | — |
| **nf-core/atacseq** | ATAC-seq chromatin accessibility | — |
| **nf-core/methylseq** | Bisulfite-seq methylation | — |
| **nf-core/taxprofiler** | Metagenomic profiling | — |
| **nf-core/viralrecon** | Viral genome assembly and variant calling | — |

### Running nf-core/sarek

nf-core/sarek is the direct production-grade equivalent of what you built manually in Modules 5 and 6 — it implements the full GATK Best Practices germline variant calling pipeline.

```bash
# Install Nextflow if not already installed
conda install -c bioconda nextflow -y

# Create a samplesheet CSV
# samplesheet.csv format:
# patient,sample,lane,fastq_1,fastq_2
cat > samplesheet.csv << 'EOF'
patient,sample,lane,fastq_1,fastq_2
PATIENT1,SAMPLE1,lane1,raw_data/SRR6821753_R1.fastq.gz,raw_data/SRR6821753_R2.fastq.gz
EOF

# Run nf-core/sarek
nextflow run nf-core/sarek \
    --input samplesheet.csv \
    --outdir sarek_results/ \
    --genome GRCh38 \
    --tools haplotypecaller \
    -profile docker \
    -resume

# For HPC without Docker, use singularity
nextflow run nf-core/sarek \
    --input samplesheet.csv \
    --outdir sarek_results/ \
    --genome GRCh38 \
    --tools haplotypecaller \
    -profile singularity \
    -resume
```

### Understanding nf-core/sarek Output

After a successful sarek run, the output directory contains:

```
sarek_results/
├── preprocessing/
│   ├── markduplicates/         ← BAM files after duplicate marking
│   └── recalibrated/           ← BQSR BAM files (analysis-ready)
├── variant_calling/
│   └── haplotypecaller/
│       └── SAMPLE1/
│           ├── SAMPLE1.vcf.gz  ← raw variant calls
│           └── SAMPLE1.vcf.gz.tbi
├── annotation/
│   └── SAMPLE1_snpeff.vcf.gz   ← annotated variants
├── reports/
│   └── multiqc/
│       └── multiqc_report.html ← aggregated QC report
└── pipeline_info/
    ├── execution_report.html   ← full run report
    └── software_versions.yml   ← exact tool versions used
```

---

## 8.6 Containerisation — Docker and Singularity

### Why Containers?

Even with conda environments, software dependency conflicts and operating system differences can cause a pipeline to behave differently on different machines. Containers solve this problem completely — a container packages the application and all its dependencies into a single portable unit that runs identically everywhere.

**Docker** — the most widely used container platform. Requires root/admin privileges to run. Common on laptops and cloud.

**Singularity** (now called Apptainer) — designed for HPC environments where users do not have root access. Can run Docker containers without modification. The standard container tool on Indian HPC clusters (PARAM Supercomputing, C-DAC systems).

### Using Docker with Nextflow

```bash
# Pull a specific tool container
docker pull broadinstitute/gatk:4.4.0.0

# Run GATK inside a container
docker run broadinstitute/gatk:4.4.0.0 \
    gatk HaplotypeCaller --version

# Nextflow handles container pulling automatically when -profile docker is used
nextflow run main.nf -profile docker
```

### Using Singularity on HPC

```bash
# Pull a Docker image as a Singularity image
singularity pull gatk.sif docker://broadinstitute/gatk:4.4.0.0

# Run a tool from a Singularity image
singularity exec gatk.sif gatk HaplotypeCaller --version

# Nextflow uses Singularity automatically
nextflow run nf-core/sarek -profile singularity
```

### Container Images for NGS Tools

Most major NGS tools have official Docker images available on:
- **BioContainers** (https://biocontainers.pro) — community-maintained containers for bioinformatics tools; used by nf-core
- **Quay.io** — hosts most BioContainers images
- **Docker Hub** — Broad Institute hosts official GATK images here

---

## 8.7 Running Pipelines on HPC and Cloud

### Running on an HPC Cluster (SLURM)

Most Indian academic institutions with bioinformatics infrastructure use SLURM (Simple Linux Utility for Resource Management) as their job scheduler. Nextflow integrates directly with SLURM — you write your pipeline once and add a SLURM profile to run it on the cluster.

```bash
# Create a SLURM submit script for Nextflow
cat > run_pipeline_slurm.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=ngs_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

module load nextflow/23.10
module load singularity/3.8

nextflow run nf-core/sarek \
    --input samplesheet.csv \
    --outdir results/ \
    --genome GRCh38 \
    --tools haplotypecaller \
    -profile singularity \
    -c nextflow_slurm.config \
    -resume
EOF

sbatch run_pipeline_slurm.sh
```

**SLURM job monitoring commands:**

```bash
squeue -u $USER          # view your running/pending jobs
squeue -j <job_id>       # view a specific job
sacct -j <job_id>        # view completed job accounting
scancel <job_id>         # cancel a job
sinfo                    # view available partitions/nodes
```

### Running on AWS (Cloud)

For large projects or when local HPC is unavailable, cloud computing is increasingly accessible. The AWS Mumbai region (`ap-south-1`) is the closest to India and offers low latency.

The key services for NGS on AWS are:
- **AWS Batch** — runs containerised jobs at scale; integrates directly with Nextflow
- **S3** — object storage for FASTQ files, reference genomes, and results
- **EC2** — virtual machines; useful for interactive analysis

```bash
# Run nf-core/sarek on AWS Batch
nextflow run nf-core/sarek \
    --input s3://my-bucket/samplesheet.csv \
    --outdir s3://my-bucket/results/ \
    --genome GRCh38 \
    --tools haplotypecaller \
    -profile awsbatch \
    -work-dir s3://my-bucket/work/ \
    -resume
```

> 💡 **Cost awareness:** Cloud computing costs money. Always estimate costs before running a large job. A 30× WGS sample through sarek on AWS costs approximately $5–15 USD depending on instance types and spot pricing. Use spot instances where possible — they are 60–80% cheaper than on-demand.

---

## 8.8 Pipeline Best Practices

### Project Directory Structure

Adopt a consistent directory structure for every project:

```
project_name/
├── README.md               ← what this project is, how to run it
├── samplesheet.csv         ← sample metadata
├── Snakefile               ← or main.nf for Nextflow
├── nextflow.config         ← configuration
├── envs/                   ← conda environment YAML files
│   └── ngs_tools.yaml
├── raw_data/               ← original FASTQ files — READ ONLY
├── reference/              ← reference genome and indices
├── results/
│   ├── qc/                 ← FastQC, MultiQC, fastp reports
│   ├── trimmed/            ← preprocessed FASTQ
│   ← aligned/             ← sorted, indexed BAM files
│   ├── variants/           ← VCF files
│   └── reports/            ← summary reports and metrics
├── logs/                   ← pipeline logs
└── scripts/                ← custom analysis scripts
```

### Version Control Everything

```bash
# Initialise a git repository for your project
git init
git add Snakefile nextflow.config envs/ scripts/ README.md samplesheet.csv
git commit -m "Initial pipeline setup"

# What NOT to commit — add to .gitignore
cat > .gitignore << 'EOF'
raw_data/
reference/
results/
logs/
work/
.nextflow/
*.bam
*.bai
*.vcf.gz
*.fastq.gz
EOF
```

### Record Software Versions

```bash
# For conda environments
conda env export > envs/environment_full.yml
conda env export --from-history > envs/environment_minimal.yml

# For Nextflow/nf-core pipelines
# software_versions.yml is generated automatically in results/pipeline_info/

# Manual version recording
echo "Tool versions used:" > logs/tool_versions.txt
fastqc --version >> logs/tool_versions.txt
fastp --version 2>&1 | head -1 >> logs/tool_versions.txt
bwa-mem2 version 2>&1 | head -1 >> logs/tool_versions.txt
samtools --version | head -1 >> logs/tool_versions.txt
gatk --version >> logs/tool_versions.txt
snpEff -version 2>&1 | head -1 >> logs/tool_versions.txt
```

### Write a Methods Section As You Go

Every parameter you choose should be documented. When writing your paper or thesis methods section, you should be able to state exactly:
- Which tool and version was used
- What parameters were used and why
- Which reference genome (including version and URL)
- Which annotation databases (including version and date downloaded)

---

## 8.9 Hands-on Exercises

### Exercise 8.1 — Install Snakemake and Run a Dry Run

```bash
conda install -c bioconda -c conda-forge snakemake -y

mkdir -p ~/ngs_workshop/module8 && cd ~/ngs_workshop/module8

# Create the Snakefile from Section 8.3
# (copy the Snakefile content from the module)
nano Snakefile

# Create required directories
mkdir -p raw_data reference results/{qc/fastqc,qc/fastp,trimmed,aligned,variants,reports,multiqc}

# Link data from previous modules
ln -s ~/ngs_workshop/module4/raw_data/SRR6821753_1.fastq.gz raw_data/SRR6821753_R1.fastq.gz
ln -s ~/ngs_workshop/module4/raw_data/SRR6821753_2.fastq.gz raw_data/SRR6821753_R2.fastq.gz
ln -s ~/ngs_workshop/module5/reference/chr20.fa reference/chr20.fa
ln -s ~/ngs_workshop/module5/reference/dbsnp138_chr20.vcf.gz reference/

# Update SAMPLES in Snakefile to just your one sample
# SAMPLES = ["SRR6821753"]

# Dry run — see what would be executed
snakemake --dry-run --cores 4

# How many steps would be run?
snakemake --dry-run --cores 4 | grep "rule" | wc -l
```

**Questions:**
1. How many steps does Snakemake plan to run?
2. In what order will the steps execute?
3. Which steps could run in parallel if you had multiple samples?

---

### Exercise 8.2 — Visualise the Pipeline DAG

```bash
# Install graphviz for visualising the DAG
conda install -c conda-forge graphviz -y

# Generate and save the DAG
snakemake --dag | dot -Tpng > pipeline_dag.png
echo "DAG saved: pipeline_dag.png"

# Also generate a rulegraph (simpler — shows rules not individual files)
snakemake --rulegraph | dot -Tpng > pipeline_rulegraph.png
echo "Rulegraph saved: pipeline_rulegraph.png"
```

Open `pipeline_rulegraph.png` to see a visual representation of your pipeline — this is excellent for presentations and methods figures.

---

### Exercise 8.3 — Run the Snakemake Pipeline

```bash
cd ~/ngs_workshop/module8

# Run the full pipeline with 4 cores
snakemake --cores 4 --use-conda

# Monitor progress — Snakemake prints each step as it starts and finishes
# If a step fails, Snakemake reports exactly which rule and sample failed

# Check outputs
ls -lh results/variants/
ls -lh results/reports/

echo "Final variant count:"
bcftools stats results/variants/SRR6821753_PASS.vcf.gz | grep "number of records"
```

---

### Exercise 8.4 — Install Nextflow and Run a Hello World

```bash
conda install -c bioconda nextflow -y
nextflow -version

# Run a simple hello world pipeline
nextflow run hello

# Run the nf-core/demo pipeline (small test dataset included)
nextflow run nf-core/demo -profile test,conda --outdir nfcore_test/
```

---

### Exercise 8.5 — Explore nf-core/sarek Documentation

Visit https://nf-co.re/sarek and answer these questions:

```
□ What variant callers does sarek support?
  (Hint: check the --tools parameter documentation)

□ What is the minimum recommended coverage for germline WGS?

□ Does sarek support somatic (tumour/normal) variant calling?

□ What containers does sarek use by default?

□ How would you run sarek on chromosome 20 only to save time?
  (Hint: look for --intervals parameter)

□ What does the sarek MultiQC report contain?
```

---

### Exercise 8.6 — Write a Methods Paragraph

Using the tools and parameters from this workshop, write a complete methods paragraph suitable for a research paper. Use the template below:

```
Raw paired-end sequencing reads were assessed for quality using FastQC
(v______) and MultiQC (v______). Adapter trimming and quality filtering
were performed using ______ (v______) with the following parameters:
______. Trimmed reads were aligned to the human reference genome
(GRCh38/hg38) using ______ (v______). Duplicate reads were marked
using ______ (v______). Base quality score recalibration was performed
using ______ with known variant sites from ______. Germline variant
calling was performed using ______ (v______) following GATK Best
Practices guidelines. Variants were filtered using hard filters with
thresholds ______. Functional annotation was performed using ______ 
(v______) with the ______ database. The full pipeline was implemented
using ______ (v______) to ensure reproducibility.
```

Fill in every blank with the tool name, version, and parameters you used in this workshop.

---

## 8.10 Where to Go From Here

You have now completed a full 2-day NGS workshop covering every step from raw FASTQ files to annotated, interpreted variant calls. Here is how to continue building on this foundation:

### Deepen Your Skills

**Bioinformatics fundamentals:**
- *Bioinformatics Data Skills* by Vince Buffalo (O'Reilly) — the best practical book for command-line bioinformatics
- *Bioinformatics Algorithms* by Compeau & Pevzner — algorithmic foundations of sequence analysis

**NGS analysis:**
- GATK Best Practices documentation — https://gatk.broadinstitute.org
- nf-core pipeline documentation — https://nf-co.re
- Bioconductor workflows — https://bioconductor.org/packages/release/BiocViews.html#___Workflow

**Statistics for genomics:**
- *Modern Statistics for Modern Biology* by Holmes & Huber — free online; excellent for R-based genomics statistics

### Practice Datasets

| Dataset | Type | Where to find |
|---|---|---|
| NA12878 (HG001) | Human WGS 30× (gold standard) | NCBI SRA, GIAB |
| 1000 Genomes Project | Population WGS | https://www.internationalgenome.org |
| TCGA | Cancer WGS/WES/RNA-seq | https://portal.gdc.cancer.gov |
| GTEx | RNA-seq across tissues | https://gtexportal.org |
| SRA public datasets | Everything | https://www.ncbi.nlm.nih.gov/sra |

### Communities and Forums

- **Biostars** (https://www.biostars.org) — bioinformatics Q&A; search before posting
- **SEQanswers** (http://seqanswers.com) — NGS-specific discussions
- **GATK forum** (https://gatk.broadinstitute.org/hc/en-us/community/topics) — official support for GATK questions
- **nf-core Slack** (https://nf-co.re/join/slack) — active community for Nextflow/nf-core help
- **Bioinformatics India** groups on LinkedIn and Twitter/X — growing community of Indian bioinformaticians

### Career Paths in Genomics and Bioinformatics in India

The skills you have learned in this workshop are in strong and growing demand across India:

- **Academic research** — NCBS, IGIB, CCMB, NIMHANS, IISc, IITs all have active genomics groups
- **Clinical genomics** — MedGenome, 4baseCare, Strand Life Sciences, Prognomics
- **Pharma and biotech** — Sun Pharma, Biocon, Serum Institute — increasingly building bioinformatics teams
- **Government initiatives** — GenomeIndia project, CSIR-IGIB, DBT-funded centres

---

> 💡 **Final takeaway:** The manual pipeline you built over these two days — FASTQ → QC → trimming → alignment → variant calling → annotation → visualisation — is exactly the same pipeline that runs in production at clinical genomics laboratories worldwide. The difference is that production pipelines run it automatically, reproducibly, and at scale using tools like Snakemake and Nextflow. You now understand every step well enough to build, debug, and interpret those production pipelines. That is a genuinely valuable and rare skill.

---

**Previous:** [Module 7 — Visualization and Interpretation](./module7_visualization.md)  
**Back to start:** [Workshop Overview](../README.md)

---

*🧬 End of NGS Workshop — 2-Day Intensive Program*  
*Thank you for participating. Happy sequencing!*
