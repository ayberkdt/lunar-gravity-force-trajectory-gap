"""Download the external data products the experiments need.

The gravity fields and SPICE kernels are public archive products and are not
redistributed here. This script fetches them into ``data/`` and verifies each
gravity file against the SHA-256 recorded in the campaign manifest, so a
mismatch is caught immediately rather than surfacing as a wrong number later.

    python fetch_data.py --list          # show what is needed, and its status
    python fetch_data.py --group lunar   # primary field only (189 MB)
    python fetch_data.py                 # everything (~360 MB)

Downloads resume-safely: a file whose digest already matches is skipped.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
MANIFEST = ROOT / "metrics" / "r16_final_experiment_manifest.json"

PDS = "https://pds-geosciences.wustl.edu"
TUDAT = ("https://raw.githubusercontent.com/tudat-team/tudat-resources"
         "/master/resource/gravity_models")
NAIF = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels"

# (group, subdirectory, filename, url, manifest key for the digest)
PRODUCTS = [
    ("lunar", "gravity_models", "jggrx_1800f_sha.tab",
     f"{PDS}/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/jggrx_1800f_sha.tab",
     "JGGRX_1800F"),
    ("lunar", "gravity_models", "gggrx_1200a_sha.tab",
     f"{PDS}/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200a_sha.tab",
     "GRGM1200A"),
    ("lunar", "gravity_models", "gggrx_1200l_sha.tab",
     f"{TUDAT}/Moon/gggrx_1200l_sha.tab", "GGGRX_1200L"),
    ("bodies", "gravity_models", "jgmess_160a_sha.tab",
     f"{TUDAT}/Mercury/jgmess_160a_sha.tab", "JGMESS_160A"),
    ("bodies", "gravity_models", "shgj180u.a01",
     f"{TUDAT}/Venus/shgj180u.a01", "SHGJ180U"),
    ("bodies", "gravity_models", "GOCO05c.txt",
     f"{TUDAT}/Earth/GOCO05c.txt", "GOCO05c"),
    ("bodies", "gravity_models", "egm96.txt",
     f"{TUDAT}/Earth/egm96.txt", "EGM96"),
    ("bodies", "gravity_models", "jgmro120d.txt",
     f"{TUDAT}/Mars/jgmro120d.txt", "JGMRO120D"),
    # SPICE kernels carry no digest in the campaign manifests; the computed
    # value is printed so it can be pinned downstream.
    ("spice", "spice_kernels", "de440s.bsp",
     f"{NAIF}/spk/planets/de440s.bsp", None),
    ("spice", "spice_kernels", "naif0012.tls",
     f"{NAIF}/lsk/naif0012.tls", None),
    ("spice", "spice_kernels", "gm_de440.tpc",
     f"{NAIF}/pck/gm_de440.tpc", None),
    ("spice", "spice_kernels", "moon_pa_de440_200625.bpc",
     f"{NAIF}/pck/moon_pa_de440_200625.bpc", None),
    ("spice", "spice_kernels", "moon_de440_250416.tf",
     f"{NAIF}/fk/satellites/moon_de440_250416.tf", None),
]

GROUPS = ("lunar", "bodies", "spice")


def recorded_digests() -> dict[str, str]:
    if not MANIFEST.exists():
        return {}
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    out = {}
    for key, meta in (manifest.get("input_products") or {}).items():
        if isinstance(meta, dict) and meta.get("sha256"):
            out[key] = meta["sha256"]
    return out


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    with urllib.request.urlopen(url) as response:
        total = response.headers.get("Content-Length")
        total = int(total) if total else 0
        done = 0
        with partial.open("wb") as handle:
            while True:
                block = response.read(1 << 20)
                if not block:
                    break
                handle.write(block)
                done += len(block)
                if total:
                    pct = 100.0 * done / total
                    print(f"\r    {done / 1048576:8.1f} MB  {pct:5.1f}%",
                          end="", flush=True)
        if total:
            print()
    shutil.move(str(partial), str(target))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", choices=GROUPS, action="append",
                        help="restrict to a group (repeatable)")
    parser.add_argument("--list", action="store_true",
                        help="report status without downloading")
    args = parser.parse_args()

    wanted = set(args.group) if args.group else set(GROUPS)
    digests = recorded_digests()
    failures = 0

    for group, subdir, name, url, key in PRODUCTS:
        if group not in wanted:
            continue
        target = DATA / subdir / name
        expected = digests.get(key) if key else None

        if target.exists():
            actual = sha256(target)
            if expected is None:
                print(f"[have]    {subdir}/{name}  sha256 {actual}")
                continue
            if actual == expected:
                print(f"[ok]      {subdir}/{name}  digest verified")
                continue
            print(f"[differs] {subdir}/{name}\n"
                  f"          expected {expected}\n"
                  f"          actual   {actual}")
            failures += 1
            continue

        if args.list:
            print(f"[missing] {subdir}/{name}  <- {url}")
            continue

        print(f"[fetch]   {subdir}/{name}")
        try:
            download(url, target)
        except Exception as exc:  # network, permissions, disk
            print(f"[error]   {name}: {exc}")
            failures += 1
            continue

        actual = sha256(target)
        if expected is None:
            print(f"[done]    {name}  sha256 {actual}")
        elif actual == expected:
            print(f"[done]    {name}  digest verified")
        else:
            print(f"[BAD]     {name}\n"
                  f"          expected {expected}\n"
                  f"          actual   {actual}")
            failures += 1

    if failures:
        print(f"\n{failures} product(s) missing or failed verification")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
