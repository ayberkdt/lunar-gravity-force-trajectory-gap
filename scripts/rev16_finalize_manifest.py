"""SHA-256 integrity manifest for the R16 transfer test.

R16 is the only campaign in the package that involves no propagation. It reads
eight distributed coefficient products and puts each through the same
field-level pipeline used for the main calibration, so the manifest indexes the
input files themselves alongside the driver, the record, and the two generated
tables.

Usage:  python rev16_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r16_final_experiment_manifest.json"

SCRIPTS = [
    "rev16_multibody_calibration.py",
    "rev16_tables.py",
    "rev16_finalize_manifest.py",
]

# r6_lp150q_transfer predates the campaign-manifest system, feeds the LP150Q
# regularization discussion, and belonged to no manifest until the reverse
# ownership scan asked; the transfer campaign is its natural home
RESULT_JSON = ["r16_multibody_calibration.json", "r6_lp150q_transfer.json"]

TABLES = ["r16_transfer_table.tex", "r16_transfer_detail_table.tex"]


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path) -> dict:
    out = {}
    for n in names:
        p = base / n
        out[n] = ({"sha256": sha(p), "bytes": p.stat().st_size}
                  if p.exists() else {"missing": True})
    return out


def index_inputs() -> dict:
    """Input coefficient products, hashed from the record the driver wrote.

    These eight files are distributed by their archives and are not shipped in
    the package, so the digest is the whole of the provenance. The driver
    records the absolute path it read on the machine that ran the campaign;
    writing that path into the manifest as the identity of the product made the
    manifest machine-bound, and the integrity check could only pass on this one
    machine. The distributed file name is the portable identity, the digest is
    the anchor, and the capture path is kept as a note about where it was read
    rather than as something a reader is expected to resolve.
    """
    src = METRICS / "r16_multibody_calibration.json"
    if not src.exists():
        return {}
    d = json.loads(src.read_text(encoding="utf-8"))
    return {f["key"]: {"body": f["body"], "center": f["center"],
                       "role": f["role"],
                       "distributed_file_name": Path(f["file"]).name,
                       "sha256": f["file_sha256"],
                       "max_degree_in_file": f["max_degree_in_file"],
                       "shipped_in_package": False,
                       "obtain_from": ("the distributing archive; verify the "
                                       "digest above after download"),
                       "capture_path_note": f["file"]}
            for f in d["fields"]}


def main() -> int:
    payload = {
        "schema": "r16_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "scope": ("R16 transfer test: the static-selection calibration repeated "
                  "on two further lunar solutions (cross-solution "
                  "reproducibility) and on five products of four other bodies "
                  "(cross-body procedural transfer). Field level only; no "
                  "trajectory is propagated and no trajectory-level claim of "
                  "the paper is extended by it."),
        "relationship_to_earlier_campaigns": (
            "additive and independent. R16 recomputes the JGGRX_1800F "
            "calibration from the archived coefficient file as an internal "
            "consistency check and reproduces the published p_spec = 2.134 and "
            "p_fit = 1.759; no earlier result is superseded."),
        "numerical_kernel": {
            "note": ("not exercised. R16 evaluates degree variances and the "
                     "tail criterion directly on the coefficient arrays; no "
                     "spherical-harmonic synthesis or integration is performed."),
        },
        "input_products": index_inputs(),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
    }
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}  manifest_sha256={payload['manifest_sha256'][:16]}")
    missing = [k for sec in ("scripts", "result_json", "generated_tables")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
