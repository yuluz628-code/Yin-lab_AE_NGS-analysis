# Yin-lab_AE_NGS-analysis

This script is designed to **automate the processing of CRISPR-Cas9 amplicon sequencing data**, including:

* **Preprocessing of raw sequencing data** (quality control, merging, and trimming)
* **Demultiplexing samples** based on barcodes
* **Alignment to reference sequences** and **mutation analysis using CRISPResso**
* **Secondary filtering and requantification** of mutation results
* **Summarization of final statistics** across all samples

This pipeline enables **fully automated analysis** from raw `FASTQ.gz` files to a comprehensive **mutation frequency summary table** with a single command.

---
