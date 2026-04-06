# el_gato Illumina baseline (Codespaces)

## Container
staphb/elgato:1.22.0

## Entry point
/usr/local/bin/el_gato.py

## Input
Paired-end Illumina FASTQs:
- ERR10879100_1.fastq.gz
- ERR10879100_2.fastq.gz

## Command
```bash
docker run --rm \
  -v "$PWD/data:/data" \
  -v "$PWD/results:/results" \
  staphb/elgato:1.22.0 \
  el_gato.py \
    --read1 /data/illumina/sample01/ERR10879100_1.fastq.gz \
    --read2 /data/illumina/sample01/ERR10879100_2.fastq.gz \
    --out /results/el_gato_illumina_baseline
``