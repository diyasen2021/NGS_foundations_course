#!/bin/bash
# =============================================================
# NGS Workshop — Automatic Environment Setup
# This runs once when a student opens the repo in Codespaces.
# They don't need to touch this file.
# =============================================================

set -e  # Stop if anything fails

echo ""
echo "======================================"
echo "  NGS Workshop Environment Setup"
echo "  This will take about 3-5 minutes..."
echo "======================================"
echo ""

# --------------------------------------------------------------
# Step 1 — Install Miniconda
# --------------------------------------------------------------
echo "[1/4] Installing Miniconda..."

wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p /root/miniconda3
rm /tmp/miniconda.sh

# Add conda to PATH permanently
echo 'export PATH="/root/miniconda3/bin:$PATH"' >> ~/.bashrc
export PATH="/root/miniconda3/bin:$PATH"

# Initialise conda
/root/miniconda3/bin/conda init bash
source ~/.bashrc

echo "    ✅ Miniconda installed"

# --------------------------------------------------------------
# Step 2 — Configure conda channels
# --------------------------------------------------------------
echo "[2/4] Configuring conda channels..."

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

echo "    ✅ Channels configured"

# --------------------------------------------------------------
# Step 3 — Install all NGS tools
# --------------------------------------------------------------
echo "[3/4] Installing NGS tools (this is the slow part)..."

conda install -y \
    fastqc \
    multiqc \
    fastp \
    trimmomatic \
    cutadapt \
    sra-tools \
    nanostat \
    nanoplot \
    python=3.10

echo "    ✅ All tools installed"

# --------------------------------------------------------------
# Step 4 — Create workshop directory structure
# --------------------------------------------------------------
echo "[4/4] Setting up workshop folders..."

mkdir -p ~/ngs_workshop/{module3,module4}/{raw_data,trimmed,fastqc_before,fastqc_after,reports,logs}

echo "    ✅ Folders created"

# --------------------------------------------------------------
# Done — print a summary of installed tools
# --------------------------------------------------------------
echo ""
echo "======================================"
echo "  Setup Complete! Tools available:"
echo "======================================"
echo ""
fastqc --version
echo "multiqc    $(multiqc --version)"
echo "fastp      $(fastp --version 2>&1 | head -1)"
echo "cutadapt   $(cutadapt --version)"
echo "trimmomatic $(trimmomatic -version 2>&1)"
echo ""
echo "Your workshop folders are ready at ~/ngs_workshop/"
echo ""
echo "To get started, open a terminal and type:"
echo "  cd ~/ngs_workshop"
echo ""
echo "Happy sequencing! 🧬"
echo "======================================"
