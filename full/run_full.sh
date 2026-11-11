#!/usr/bin/env bash
# Whole-genome build: split by chromosome, build each, merge.
# Assumes full/raw/{gff3.gz,genome.fa.gz,pep.fa.gz} are downloaded.
set -euo pipefail
cd "$(dirname "$0")/.."

python3 full/split.py
python3 full/build_all.py "$@"
python3 full/merge.py
