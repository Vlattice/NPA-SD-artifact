"""Download and prepare inputs: VCFs, the pedigree, and per-region extracts.

Large VCFs are cached and downloaded with ``wget -c`` so an interrupted transfer
resumes instead of restarting. Downloads are verified against the server's
Content-Length before a file is treated as complete.
"""

from __future__ import annotations

import os
import subprocess
import time
import urllib.request
from typing import Optional

import pandas as pd

from . import config


def remote_size(url: str) -> Optional[int]:
    """Content-Length of a remote file, or None if the server does not report it."""
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length else None
    except Exception:
        return None


def download(url: str, dest: str, expected: Optional[int] = None,
             max_tries: int = 10, quiet: bool = False) -> bool:
    """Resumable download. Returns True once the file is complete on disk."""
    if expected is None:
        expected = remote_size(url)
    for attempt in range(1, max_tries + 1):
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if expected and have == expected:
            return True
        if not quiet:
            print(f"  {os.path.basename(dest)}: {have:,} B"
                  + (f" / {expected:,}" if expected else ""))
        subprocess.run(["wget", "-c", "-q", "-T", "60", "--tries", "3",
                        "-O", dest, url], check=False)
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if expected and have == expected:
            return True
        if expected is None and have > 0:
            return True
        time.sleep(min(30, 2 ** attempt))
    return False


def fetch_vcf(chrom: str, cache_dir: str) -> Optional[str]:
    """Download a chromosome VCF and its index into ``cache_dir``; return the path."""
    os.makedirs(cache_dir, exist_ok=True)
    name = config.vcf_filename(chrom)
    url = config.VCF_BASE_URL + name
    dest = os.path.join(cache_dir, name)
    if not download(url, dest):
        return None
    if not os.path.exists(dest + ".tbi"):
        if not download(url + ".tbi", dest + ".tbi"):
            subprocess.run(["bcftools", "index", "-t", dest], check=False)
    return dest


def extract_region(vcf_path: str, chrom: str, start: int, end: int,
                   out_path: str) -> str:
    """Extract one region with bcftools and index the result."""
    subprocess.run(
        f"bcftools view -r {chrom}:{start}-{end} -Oz -o {out_path} {vcf_path} "
        f"&& bcftools index -t {out_path}",
        shell=True, check=True,
    )
    return out_path


def load_trios(pedigree_path: str) -> list[dict]:
    """Parse the 1000 Genomes pedigree into complete mother-father-child trios.

    Column names vary between pedigree releases, so the sample, father, mother and
    population columns are located by name rather than by position.
    """
    ped = pd.read_csv(pedigree_path, sep=r"\s+", engine="python", dtype=str)

    def find(*needles):
        for col in ped.columns:
            low = str(col).lower()
            if any(n in low for n in needles):
                return col
        return None

    c_sample = find("sample", "individual") or ped.columns[0]
    c_father = find("father", "paternal")
    c_mother = find("mother", "maternal")
    c_pop = find("population", "pop")
    if not (c_father and c_mother):
        raise ValueError(f"could not locate parent columns in {list(ped.columns)}")

    def missing(value):
        return str(value).strip() in {"0", "", ".", "nan", "NA", "None"}

    trios = []
    for _, row in ped.iterrows():
        father, mother = row[c_father], row[c_mother]
        if missing(father) or missing(mother):
            continue
        trios.append({
            "child": str(row[c_sample]),
            "father": str(father),
            "mother": str(mother),
            "pop": str(row[c_pop]) if c_pop else "",
        })
    return trios
