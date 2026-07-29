"""Trimmed package initializer for the paper source snapshot.

The upstream Lunaris initializer eagerly imports the full propagation stack,
which is outside the scope of this snapshot. Only the submodules the published
experiments import are carried here; the modules themselves are verbatim.
"""

from __future__ import annotations

__all__: list[str] = []
