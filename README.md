# mSPOT

> **Note on Data Availability**  
> Raw `pod5` signal files are too large for GitHub storage. We provide a short FASTA file (`calls.fa`) for testing the downstream bioinformatics pipeline. Reviewers can execute **Part B, Steps 4–12** directly to verify data processing.



## Overview

This repository contains the computational design scripts and bioinformatics pipeline for **mSPOT**, a multiplexed DNAzyme sensing platform coupled with Oxford Nanopore Technologies (ONT) sequencing for parallel readout.

- **Part A** – Python scripts for designing multiplexed DNAzyme catalytic domains and nanopore barcode sequences with constrained sequence properties (GC content, distance metrics, homopolymer limits).
- **Part B** – Shell commands for basecalling, demultiplexing, alignment, and read-length analysis of nanopore sequencing data.




## Part A: Sequence Design

### 1. DNAzyme Domain Design

**File:** `script/DomainsDesign.py`

Generates DNAzyme catalytic domains with constrained sequence length, GC content, and **minimum Hamming distance** to minimize cross-reactivity among multiplexed DNAzymes.

**Parameters used in this study:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `sequence_length` | 6 nt | Domain length |
| `symbols` | A, G, C, T | Nucleotide alphabet |
| `min_distance` | 4 | Minimum Hamming distance |
| `gc_min` / `gc_max` | 2 / 4 | Total G+C count per domain |

**Usage:**

```bash
python script/DomainsDesign.py
```

**Output:** 56 sequences satisfying all constraints (Hamming distance ≥ 4, GC count 2–4).

---

### 2. Nanopore Barcode Design

**File:** `script/BarcodeDesign.py`

Generates 24-nt barcode sequences for nanopore sequencing with constraints on **Levenshtein distance**, **GC content**, and **homopolymer length** to ensure robust barcode discrimination and basecalling accuracy.

**Parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `num_sequences` | 112 | Number of barcodes to generate |
| `seq_length` | 24 nt | Barcode length |
| `min_levenshtein_distance` | 10 | Minimum Levenshtein distance between any two barcodes |
| `max_homopolymer` | 2 | Maximum allowed homopolymer run (e.g., `AA` allowed, `AAA` forbidden) |
| `gc_content_range` | 0.4 – 0.6 | Fractional GC content |

**Usage:**

```bash
python script/BarcodeDesign.py
```



---

## Part B: Nanopore Data Processing

### Prerequisites

| Tool | Version | Installation |
|------|---------|--------------|
| [Dorado](https://github.com/nanoporetech/dorado) | ≥ 0.7.0 | `conda install -c nanoporetech dorado` or pre-built binary |
| [SeqKit](https://github.com/shenwei356/seqkit) | 2.6.1 | `conda install -c bioconda seqkit` |
| [BWA](https://github.com/lh3/bwa) | 0.7.17-r1188 | `conda install -c bioconda bwa` |
| [Samtools](https://github.com/samtools/samtools) | 1.13 | `conda install -c bioconda samtools` |

All tools are available via Bioconda and install in minutes.

---

### Basecalling & Demultiplexing

```bash
# 1. Basecalling (SUP model, all CUDA devices, native barcode kit, no trimming)
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/ \
    --kit-name SQK-NBD114-24 \
    --device cuda:all \
    --no-trim > calls_notrim.bam

# 2. Demultiplexing & FASTQ generation
dorado demux -t 200 \
    --emit-summary \
    --emit-fastq \
    --kit-name SQK-NBD114-24 \
    -v \
    --output-dir calls.fq \
    calls_notrim.bam

# 3. Convert FASTQ to FASTA
seqkit fq2fa calls.fq > calls.fa
```

---

### Read Length Filtering

```bash
# 4. Filter reads by length (55–77 nt)
seqkit seq -m 55 -M 77 calls.fa > calls_m55M77.fa
```

---

### Mapping & Counting

```bash
# 5. Align filtered reads to reference (ONT-optimized preset)
bwa mem -k 10 -x ont2d ref/ref.fa calls_m55M77.fa -t 200 > calls_m55M77.sam

# 6. Sort and convert to BAM
samtools sort -@ 200 -O bam -o calls_m55M77.sorted.bam calls_m55M77.sam

# 7. Index BAM
samtools index calls_m55M77.sorted.bam

# 8. Generate per-reference read counts
samtools idxstats calls_m55M77.sorted.bam > calls_m55M77.txt
```

---

### Read Length Distribution

```bash
# 9. Extract reads ≤ 150 nt for length distribution analysis
seqkit seq -M 150 calls.fa > calls_M150.fa

# 10. Extract lengths to tabular format
seqkit fx2tab -j 30 -l -n -i -H calls_M150.fa | cut -f 2 > Length_calls_M150.txt

# 11. Compute length-frequency table
awk 'NR > 1 { print $1 }' Length_calls_M150.txt | \
    sort | uniq -c | \
    awk '{ print $2, $1 }' > Length_distribution_calls_M150.txt

# 12. Import Length_distribution_calls_M150.txt into Excel
#     to calculate percentage for each read length.
```

---

## Software Versions

| Software | Version | URL |
|----------|---------|-----|
| SeqKit | 2.6.1 | https://github.com/shenwei356/seqkit |
| BWA | 0.7.17-r1188 | https://github.com/lh3/bwa |
| Samtools | 1.13 | https://github.com/samtools/samtools |

Installation guides are available via their respective GitHub repositories. All tools install rapidly and are ready for immediate use after installation.

---

## Citation

If you use the mSPOT pipeline or design scripts in your research, please cite:

> *[Manuscript citation to be added upon publication]*

---

