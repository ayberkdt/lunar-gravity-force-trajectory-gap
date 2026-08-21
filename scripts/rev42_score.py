"""R42 scoring: R37's scorer, pointed at R42's record.

The scoring convention is not restated here and not reimplemented here. This
file imports rev37_score, rebinds the two paths it reads and writes, and calls
its main. Any difference between the R37 verdict and the R42 verdict is
therefore a difference in the rows, never in the rule -- including the scoring
amendment (sign agreement among resolved comparisons only), which is applied by
the imported code exactly as it was in R37.

rev37_score.py is sealed under the R37 manifest and is imported, never edited.
Its own record and verdict files are untouched: this writes
metrics/r42_panel_verdict.json.

Usage:  python rev42_score.py [--json]
"""

from __future__ import annotations

import json
from pathlib import Path

import rev37_score as scorer

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

scorer.RECORD = METRICS / "r42_variational_completion.json"
scorer.OUT = METRICS / "r42_panel_verdict.json"


def main() -> int:
    rc = scorer.main()
    if rc == 0 and scorer.OUT.exists():
        out = json.loads(scorer.OUT.read_text(encoding="utf-8"))
        out["schema"] = "r42_panel_verdict_v1"
        out["scored_by"] = ("rev37_score.py, imported verbatim and pointed at "
                            "the R42 record; the scoring rule is R37's")
        scorer.OUT.write_text(json.dumps(out, indent=2) + "\n",
                              encoding="utf-8")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
