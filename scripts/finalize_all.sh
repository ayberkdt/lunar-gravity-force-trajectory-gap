#!/bin/sh
# Turn whatever the queue managed to propagate into readable results.
#
# This exists so the campaign does not depend on anyone being awake when it
# ends. Every step reads the records on disk and reports what is there: a
# population that never ran is reported as not run, a budget that never
# finished is not indexed, and nothing here invents a number for a stage that
# did not happen.
#
# Safe to run more than once, and safe to run early.
cd "$(dirname "$0")" || exit 1

echo "=== finalize at $(date) ==="

echo "--- design C verdict ---"
python rev29_verdict.py
echo "--- design C tables ---"
python rev29_tables.py

echo "--- operational elliptical (outside the sampled box) ---"
python rev30_verdict.py --registry r31
python rev30_tables.py --registry r31

echo "--- geometry strata ---"
python rev30_verdict.py
python rev30_tables.py

echo "--- manifests ---"
python rev29_finalize_manifest.py
python rev30_finalize_manifest.py --registry r31
python rev30_finalize_manifest.py

echo "--- integrity gate ---"
python check_manifest_integrity.py

echo "=== finalize done at $(date) ==="
