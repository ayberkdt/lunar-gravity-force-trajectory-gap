from __future__ import annotations

import unittest

from tools.check_frozen_changes import Change, classify, parse_name_status


class ParseNameStatusTests(unittest.TestCase):
    def test_parses_regular_and_rename_rows(self) -> None:
        changes = parse_name_status(
            "A\tdocs/new.md\nM\tscripts/old.py\nR100\ta.py\tb.py\n"
        )
        self.assertEqual(
            changes,
            [
                Change("A", ("docs/new.md",)),
                Change("M", ("scripts/old.py",)),
                Change("R100", ("a.py", "b.py")),
            ],
        )

    def test_rejects_malformed_rows(self) -> None:
        with self.assertRaises(ValueError):
            parse_name_status("R100\tonly-one-path")


class ArchiveBoundaryTests(unittest.TestCase):
    def assert_allowed(self, change: Change) -> None:
        _allowed, blocked = classify([change])
        self.assertEqual(blocked, [])

    def assert_blocked(self, change: Change) -> None:
        _allowed, blocked = classify([change])
        self.assertEqual(len(blocked), 1)

    def test_allows_new_documentation_and_tools(self) -> None:
        self.assert_allowed(Change("A", ("docs/ARCHITECTURE.md",)))
        self.assert_allowed(Change("A", ("tools/new_check.py",)))

    def test_allows_regenerating_maintained_docs(self) -> None:
        self.assert_allowed(Change("M", ("environments/README.md",)))
        self.assert_blocked(
            Change("M", ("environments/field-and-timing-py312.txt",))
        )

    def test_allows_new_campaign_artifact(self) -> None:
        self.assert_allowed(Change("A", ("metrics/r52_result.json",)))
        self.assert_allowed(Change("A", ("scripts/rev52_campaign.py",)))

    def test_blocks_existing_evidence_edit_or_delete(self) -> None:
        self.assert_blocked(Change("M", ("metrics/r14_preregistration.json",)))
        self.assert_blocked(Change("D", ("scripts/rev14_budget_pareto.py",)))

    def test_blocks_rename_touching_evidence(self) -> None:
        self.assert_blocked(
            Change("R100", ("scripts/rev14_budget_pareto.py",
                            "scripts/budget_pareto.py"))
        )

    def test_blocks_every_source_snapshot_change(self) -> None:
        self.assert_blocked(
            Change("M", ("src/lunaris/physics/spherical_harmonics.py",))
        )
        self.assert_blocked(Change("A", ("src/lunaris/new_module.py",)))

    def test_blocks_verification_contract_change(self) -> None:
        self.assert_blocked(Change("M", (".gitattributes",)))
        self.assert_blocked(Change("M", ("verify_snapshot.py",)))
        self.assert_blocked(Change("M", ("docs/SOURCE_MANIFEST.md",)))


if __name__ == "__main__":
    unittest.main()
