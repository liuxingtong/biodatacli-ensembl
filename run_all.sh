#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

python3 fetch_data.py
python3 make_enrichment.py
python3 crosswalk/build.py
python3 -m ensembl_adapter.build \
  --gff3 data/raw/Homo_sapiens.GRCh38.116.chromosome.20.gff3.gz \
  --genome data/raw/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz \
  --pep data/raw/Homo_sapiens.GRCh38.pep.all.fa.gz \
  --biomart data/raw/biomart \
  --enrichment data/raw/xref_enrichment.json \
  --capabilities capabilities.json \
  --db data/index/ensembl_chr20_r116.sqlite \
  --stats data/index/build_stats.json
python3 prototype_4dnk.py
python3 verify.py
python3 acceptance_m1.py
python3 contract/conformance.py
python3 bench/run.py
python3 bridge/smoke.py
python3 -m unittest discover -s tests -v
