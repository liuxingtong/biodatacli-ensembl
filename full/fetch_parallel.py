#!/usr/bin/env python3
"""Parallel range downloader for the whole-genome inputs.

ftp.ensembl.org throttles a single TCP stream to ~15 KB/s from this network, but
the cap is per connection: eight streams reach ~120 KB/s.  This downloads each
file in fixed-size chunks across a small thread pool, resumes by skipping
complete chunks, and concatenates when done.

    python3 full/fetch_parallel.py            # all inputs
    python3 full/fetch_parallel.py --workers 16
"""

from __future__ import annotations

import argparse
import os
import sys
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "full", "raw")

FTP = "https://ftp.ensembl.org/pub/release-116"
FILES = [
    (f"{FTP}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz", "pep.fa.gz", 23_319_936),
    (f"{FTP}/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz", "gff3.gz", 108_378_077),
    (f"{FTP}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
     "genome.fa.gz", 881_964_081),
]

CHUNK = 256 * 1024
_lock = threading.Lock()
_progress = {"bytes": 0, "started": time.time()}


def remote_size(url: str) -> int:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        return int(response.headers["Content-Length"])


def fetch_chunk(url: str, path: str, start: int, end: int, retries: int = 5) -> int:
    expected = end - start + 1
    if os.path.exists(path) and os.path.getsize(path) == expected:
        return expected
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(request, timeout=120) as response:
                data = response.read()
            if len(data) != expected:
                raise IOError(f"short read {len(data)}/{expected}")
            temporary = path + ".tmp"
            with open(temporary, "wb") as handle:
                handle.write(data)
            os.replace(temporary, path)
            return expected
        except Exception as exc:  # noqa: BLE001 - retry network errors
            if attempt == retries:
                raise
            time.sleep(1.5 * attempt)
    return 0


def download(url: str, name: str, expected_size: int, workers: int) -> str:
    destination = os.path.join(RAW, name)
    if os.path.exists(destination) and os.path.getsize(destination) == expected_size:
        print(f"[skip] {name} already complete")
        return destination
    size = remote_size(url) or expected_size
    chunk_dir = destination + ".chunks"
    os.makedirs(chunk_dir, exist_ok=True)
    ranges = [(index * CHUNK, min((index + 1) * CHUNK - 1, size - 1))
              for index in range((size + CHUNK - 1) // CHUNK)]
    print(f"[get ] {name}: {size:,} bytes in {len(ranges)} chunks, {workers} workers", flush=True)

    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(fetch_chunk, url, os.path.join(chunk_dir, f"{index:05d}.part"), start, end): index
            for index, (start, end) in enumerate(ranges)
        }
        for future in as_completed(futures):
            done += 1
            if done % 20 == 0 or done == len(ranges):
                elapsed = time.time() - _progress["started"]
                rate = done * CHUNK / max(1, elapsed) / 1024
                print(f"    {name}: {done}/{len(ranges)} chunks  ~{rate:.0f} KB/s", flush=True)
            future.result()

    temporary = destination + ".assembling"
    with open(temporary, "wb") as out:
        for index in range(len(ranges)):
            part = os.path.join(chunk_dir, f"{index:05d}.part")
            with open(part, "rb") as handle:
                while True:
                    block = handle.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
    actual = os.path.getsize(temporary)
    if actual != size:
        raise SystemExit(f"{name}: assembled {actual} != {size}")
    os.replace(temporary, destination)
    for index in range(len(ranges)):
        os.remove(os.path.join(chunk_dir, f"{index:05d}.part"))
    os.rmdir(chunk_dir)
    print(f"[done] {name}: {actual:,} bytes", flush=True)
    return destination


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args(argv)
    os.makedirs(RAW, exist_ok=True)
    for url, name, size in FILES:
        download(url, name, size, args.workers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
