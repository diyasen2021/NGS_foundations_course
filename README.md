# 🧬 NGS Workshop — Next Generation Sequencing
### A 2-Day Intensive Program for UG & PG Students in India

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/YOUR_GITHUB_USERNAME/YOUR_REPO_NAME)

> **Replace `YOUR_GITHUB_USERNAME` and `YOUR_REPO_NAME` in the badge link above with your actual GitHub username and repository name before sharing with students.**

---

## 🚀 Quick Start for Students

**Step 1** — Click the green button above **"Open in GitHub Codespaces"**  
**Step 2** — Wait ~3-5 minutes while the environment sets up automatically  
**Step 3** — Open a terminal and type `cd ~/ngs_workshop`  
**Step 4** — You're ready! All tools are pre-installed ✅

> You will need a free GitHub account. Sign up at [github.com](https://github.com) if you don't have one.

---

## 📋 Workshop Overview

This workshop covers the full NGS analysis pipeline from raw sequencing data to biological interpretation, using real datasets and industry-standard tools.

| | Day 1 — Fundamentals & Raw Data Analysis |
|---|---|
| Module 1 | [Introduction to Genomics and Sequencing Technologies](Lessons/module1_introduction_to_genomics.md) |
| Module 2 | [NGS Workflow and Data Generation](Lessons/module2_ngs_workflow.md) |
| Module 3 | [Quality Control of Raw Sequencing Data 🧪](Lessons/module3_quality_control.md) |
| Module 4 | [Data Preprocessing and Cleaning 🧪](Lessons/module4_preprocessing.md) |

| | Day 2 — Advanced NGS Data Analysis |
|---|---|
| Module 5 | [Sequence Alignment and Mapping 🧪](Lessons/module5_alignment.md) |
| Module 6 | [Variant Detection and Analysis 🧪](Lessons/module6_variant_calling.md) |
| Module 7 | [Visualization and Interpretation](Lessons/module7_visualization.md) |
| Module 8 | [Advanced NGS Pipelines and Workflow Management](Lessons/module8_pipelines.md) |

🧪 = Hands-on practical session

---

## 🛠️ Tools Covered

| Tool | Purpose | Module |
|---|---|---|
| FastQC | Raw data quality control | 3 |
| MultiQC | Aggregate QC reports | 3 |
| fastp | Read trimming and preprocessing | 4 |
| Trimmomatic | Read trimming (alternative) | 4 |
| Cutadapt | Adapter trimming | 4 |
| BWA-MEM2 | DNA sequence alignment | 5 |
| HISAT2 | RNA-seq alignment | 5 |
| SAMtools | BAM file processing | 5 |
| GATK | Variant calling | 6 |
| IGV | Genome visualisation | 7 |
| Nextflow / nf-core | Pipeline management | 8 |

---

## 💻 System Requirements

You only need a **web browser** and a **GitHub account** — everything else runs in the cloud via GitHub Codespaces.

| Requirement | Details |
|---|---|
| Browser | Chrome, Firefox, Edge, or Safari |
| GitHub account | Free — [sign up here](https://github.com/signup) |
| Internet connection | Required throughout the workshop |
| Local software | None needed |

---

## 🗂️ Repository Structure

```
project/
├── .devcontainer/
│   ├── devcontainer.json     ← Codespaces configuration
│   └── setup.sh              ← Auto-installs all tools on startup
├── Lessons/
│   ├── module1_introduction_to_genomics.md
│   ├── module2_ngs_workflow.md
│   ├── module3_quality_control.md
│   ├── module4_preprocessing.md
│   ├── module5_alignment.md
│   ├── module6_variant_calling.md
│   ├── module7_visualization.md
│   └── module8_pipelines.md
└── README.md
```

---

## ❓ Frequently Asked Questions

**Q: Do I need to install anything on my laptop?**  
No. Everything runs in GitHub Codespaces in your browser. Windows, Mac, and Linux laptops all work equally well.

**Q: Will my work be saved?**  
Yes — your Codespace persists between sessions as long as you don't delete it. Files you create inside `~/ngs_workshop/` will be there when you come back.

**Q: My Codespace timed out — is my work lost?**  
No. Codespaces automatically pause after 30 minutes of inactivity but your files are saved. Just reopen the Codespace from [github.com/codespaces](https://github.com/codespaces) and everything will be exactly as you left it.

**Q: How much does Codespaces cost?**  
Free GitHub accounts include 60 hours of Codespaces per month — more than enough for this 2-day workshop.

**Q: I want to keep practising at home after the workshop — what do I use?**  
Your Codespace stays active after the workshop. Alternatively, Windows users can install WSL2 (Windows Subsystem for Linux) to run the same tools locally — ask your instructor for the setup guide.

---

## 👩‍🏫 For Instructors

### Setting up the Codespaces button

After creating your GitHub repo, update the badge link at the top of this README:

```
Replace:  YOUR_GITHUB_USERNAME  →  your actual GitHub username
Replace:  YOUR_REPO_NAME        →  your actual repository name
```

For example, if your GitHub username is `drpriya` and your repo is `ngs-workshop-2025`:
```markdown
[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/drpriya/ngs-workshop-2025)
```

### Sharing with students

Send students this single link — clicking it takes them directly to open the repo in Codespaces:
```
https://codespaces.new/YOUR_GITHUB_USERNAME/YOUR_REPO_NAME
```

### Pre-loading the test dataset

To save time on the day, you can add a data download step to `setup.sh`. Add the following lines at the end of the file:

```bash
# Download test dataset (500,000 reads — ~200 MB)
cd ~/ngs_workshop/module3/raw_data
fastq-dump --split-files --gzip -X 500000 SRR6821753
```

This means the dataset is ready and waiting when students arrive at Module 3.

---

## 📚 Further Reading

- [GATK Best Practices Documentation](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)
- [nf-core Pipelines](https://nf-co.re/pipelines)
- [Bioconductor Workflows](https://bioconductor.org/packages/release/BiocViews.html#___Workflow)
- [EMBL-EBI Training](https://www.ebi.ac.uk/training/)
- [Biostars Community Forum](https://www.biostars.org/)

---

*Workshop materials developed for UG & PG students. All modules are written in Markdown and can be read directly on GitHub or downloaded for offline use.*
