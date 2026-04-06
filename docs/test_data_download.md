# Downloading test FASTQ data

This repository does not store raw sequencing data in Git.
Instead, test datasets are downloaded on demand using public
SRA accessions.

All downloads are performed using the official **NCBI SRA Toolkit
Docker container**, so no local installation is required.

---

## Test sample accessions

The initial test sample has been sequenced on both platforms:

| Sample ID | Platform  | Accession     |
|-----------|-----------|---------------|
| sample01  | Illumina  | ERR10879100   |
| sample01  | ONT       | ERR10879072   |

These mappings are also recorded in:

## Download script

A helper script is provided to download the test data
into the expected directory structure.

From the repository root:

```bash
bash scripts/download_test_data.sh


data/
├── illumina/
│   └── sample01/
│       ├── ERR10879100_1.fastq.gz
│       └── ERR10879100_2.fastq.gz
└── ont/
    └── sample01/
        └── ERR10879072.fastq.gz
