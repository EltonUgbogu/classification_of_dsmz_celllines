#!/usr/bin/env python3
"""
Build a breast cancer cell-line subtype reference table.

Only cell lines present in the DSMZ CellDive-style annotation are retained.
External subtype sources are used only when they overlap with this DSMZ set.

Outputs:
    subtype_reference_outputs/breast_subtype_sources.long.tsv
    subtype_reference_outputs/breast_subtype_sources.wide.tsv
    subtype_reference_outputs/breast_subtype_sources.consensus.tsv
    subtype_reference_outputs/breast_subtype_conflicts.tsv
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Optional

import pandas as pd


OUTPUT_DIR = Path("subtype_reference_outputs")


# ---------------------------------------------------------------------
# 1. Cell-line name normalisation
# ---------------------------------------------------------------------

NAME_ALIASES = {
    "SKBR3": "SK-BR-3",
    "SKBR-3": "SK-BR-3",
    "BT474": "BT-474",
    "T47D": "T-47D",
    "MCF7": "MCF-7",
    "MDAMB231": "MDA-MB-231",
    "MDAMB453": "MDA-MB-453",
    "MDAMB468": "MDA-MB-468",
    "HS578T": "HS-578T",
    "CAL51": "CAL-51",
    "CAL148": "CAL-148",
    "EFM19": "EFM-19",
    "EFM192A": "EFM-192A",
    "HCC1143": "HCC-1143",
    "HCC1599": "HCC-1599",
    "HCC1937": "HCC-1937",
}


def normalise_cell_line(name: str) -> str:
    """Return a consistent cell-line name."""
    cleaned = name.strip().replace(" ", "").replace("_", "-")
    key = cleaned.upper().replace("-", "")
    return NAME_ALIASES.get(key, cleaned)


# ---------------------------------------------------------------------
# 2. DSMZ CellDive-style reference annotation
# ---------------------------------------------------------------------

DSMZ_CELLDIVE = {
    "CAL-120": "B",
    "CAL-148": "L",
    "CAL-51": "B",
    "CAL-85-1": "B",
    "COLO-824": "B",
    "DU-4475": "L",
    "EFM-19": "L",
    "EFM-192A": "L",
    "EFM-192B": "L",
    "EFM-192C": "L",
    "ETCC-006": "B",
    "ETCC-007": "L",
    "EVSA-T": "L",
    "HCC-1143": "B",
    "HCC-1599": "B",
    "HCC-1937": "B",
    "HDQ-P1": "L",
    "HS-578T": "L",
    "IPH-926": "L",
    "JIMT-1": "L",
    "KPL-1": "L",
    "MCF-7": "L",
    "MDA-MB-231": "B",
    "MDA-MB-453": "L",
    "MDA-MB-468": "B",
    "MFM-223": "L",
    "SK-BR-3": "L",
    "T-47D": "L",
    "BT-474": "L",
}

DSMZ_CELL_LINES = {normalise_cell_line(cell_line) for cell_line in DSMZ_CELLDIVE}


# ---------------------------------------------------------------------
# 3. External source-specific subtype annotations
# ---------------------------------------------------------------------

IHC_MOLECULAR_CLASSIFICATION = {
    "MCF-7": "Luminal A",
    "BT-474": "Luminal B",
    "SK-BR-3": "HER2 over-expression",
    "MDA-MB-231": "Basal",
    "MDA-MB-468": "Basal",
    "HS-578T": "Basal",
    "MDA-MB-453": "Unclassified",
}

LUMINAL_BASAL_A_B_SOURCE = {
    "BT-474": "L",
    "CAL-51": "B",
    "EFM-19": "L",
    "EFM-192A": "L",
    "HCC-1143": "A",
    "HCC-1599": "A",
    "HCC-1937": "A",
    "HS-578T": "B",
    "MCF-7": "L",
    "MDA-MB-231": "B",
    "MDA-MB-453": "L",
    "MDA-MB-468": "A",
    "SK-BR-3": "L",
    "T-47D": "L",
}

RECEPTOR_TNBC_SOURCE = {
    "BT-474": "Luminal HER2",
    "MDA-MB-453": "HER2+",
    "SK-BR-3": "HER2+",
    "T-47D": "Luminal A",
    "MCF-7": "Luminal A",
    "CAL-148": "TNBCA",
    "MDA-MB-231": "TNBCB",
    "HS-578T": "TNBCB",
}


SOURCES = {
    "DSMZ_CellDive_B_L": DSMZ_CELLDIVE,
    "IHC_molecular_classification": IHC_MOLECULAR_CLASSIFICATION,
    "Luminal_Basal_A_B_source": LUMINAL_BASAL_A_B_SOURCE,
    "Receptor_TNBC_source": RECEPTOR_TNBC_SOURCE,
}


# ---------------------------------------------------------------------
# 4. Harmonisation to simplified B/L labels
# ---------------------------------------------------------------------

def collapse_to_bl(label: str) -> Optional[str]:
    """
    Collapse source-specific subtype labels to the simplified B/L scheme.

    B = basal-like / triple-negative-related
    L = luminal / HER2-related in the simplified DSMZ-style comparison
    None = unclassified, missing, or unsuitable for B/L comparison
    """
    value = label.strip().lower()

    if value in {"b", "basal", "basal-like", "tnbca", "tnbcb", "tnbc"}:
        return "B"

    if value in {
        "l",
        "luminal",
        "luminal a",
        "luminal b",
        "luminal her2",
        "her2+",
        "her2 over-expression",
        "her2 overexpression",
    }:
        return "L"

    if value == "a":
        # In luminal/basal A/B sources, A means basal A.
        return "B"

    if value in {"na", "n/a", "unclassified", "mel"}:
        return None

    raise ValueError(f"Unrecognised subtype label: {label!r}")


def subtype_family(label: str) -> str:
    """Return a descriptive subtype family."""
    value = label.strip().lower()

    if value in {"b", "basal", "basal-like", "a"}:
        return "basal-like"

    if value in {"tnbca", "tnbcb", "tnbc"}:
        return "triple-negative-related"

    if value in {"l", "luminal", "luminal a", "luminal b"}:
        return "luminal"

    if value in {"luminal her2", "her2+", "her2 over-expression", "her2 overexpression"}:
        return "HER2-related"

    if value in {"na", "n/a", "unclassified"}:
        return "unclassified"

    if value == "mel":
        return "melanoma-like_or_non-breast-like"

    return "other"


# ---------------------------------------------------------------------
# 5. Table construction
# ---------------------------------------------------------------------

def build_long_table(sources: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    rows = []

    for source_name, annotations in sources.items():
        for raw_cell_line, source_label in annotations.items():
            cell_line = normalise_cell_line(raw_cell_line)

            # Keep only cell lines present in the DSMZ CellDive-style annotation.
            if cell_line not in DSMZ_CELL_LINES:
                continue

            rows.append(
                {
                    "cell_line": cell_line,
                    "source": source_name,
                    "source_subtype": source_label,
                    "subtype_family": subtype_family(source_label),
                    "collapsed_B_L": collapse_to_bl(source_label),
                    "present_in_DSMZ_CellDive_style_annotation": True,
                }
            )

    return (
        pd.DataFrame(rows)
        .sort_values(["cell_line", "source"])
        .reset_index(drop=True)
    )


def build_wide_table(long_df: pd.DataFrame) -> pd.DataFrame:
    source_wide = long_df.pivot_table(
        index="cell_line",
        columns="source",
        values="source_subtype",
        aggfunc="first",
    )

    collapsed_wide = long_df.pivot_table(
        index="cell_line",
        columns="source",
        values="collapsed_B_L",
        aggfunc="first",
    ).add_suffix("__collapsed_B_L")

    wide_df = pd.concat([source_wide, collapsed_wide], axis=1)
    return wide_df.reset_index().sort_values("cell_line")


def majority_label(labels: Iterable[Optional[str]]) -> Optional[str]:
    usable = [label for label in labels if pd.notna(label) and label is not None]

    if not usable:
        return None

    counts = Counter(usable)
    most_common = counts.most_common()

    if len(most_common) > 1 and most_common[0][1] == most_common[1][1]:
        return "ambiguous"

    return most_common[0][0]


def build_consensus_table(long_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for cell_line, group in long_df.groupby("cell_line"):
        labels = group["collapsed_B_L"].tolist()
        observed_labels = sorted({label for label in labels if pd.notna(label)})

        rows.append(
            {
                "cell_line": cell_line,
                "dsmz_celldive_label": DSMZ_CELLDIVE[cell_line],
                "n_sources": group["source"].nunique(),
                "source_labels": "; ".join(
                    f"{row.source}={row.source_subtype}"
                    for row in group.itertuples(index=False)
                ),
                "collapsed_labels_observed": ",".join(observed_labels),
                "consensus_B_L": majority_label(labels),
                "has_B_L_conflict": len(observed_labels) > 1,
            }
        )

    return (
        pd.DataFrame(rows)
        .sort_values(["has_B_L_conflict", "cell_line"], ascending=[False, True])
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------
# 6. Main execution
# ---------------------------------------------------------------------

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    long_df = build_long_table(SOURCES)
    wide_df = build_wide_table(long_df)
    consensus_df = build_consensus_table(long_df)
    conflict_df = consensus_df[consensus_df["has_B_L_conflict"]].copy()

    long_df.to_csv(OUTPUT_DIR / "breast_subtype_sources.long.tsv", sep="\t", index=False)
    wide_df.to_csv(OUTPUT_DIR / "breast_subtype_sources.wide.tsv", sep="\t", index=False)
    consensus_df.to_csv(OUTPUT_DIR / "breast_subtype_sources.consensus.tsv", sep="\t", index=False)
    conflict_df.to_csv(OUTPUT_DIR / "breast_subtype_conflicts.tsv", sep="\t", index=False)

    print("Retained DSMZ CellDive-style cell lines:", len(DSMZ_CELL_LINES))
    print("Rows in long source table:", len(long_df))
    print("\nCollapsed B/L consensus summary:")
    print(consensus_df["consensus_B_L"].value_counts(dropna=False).to_string())

    if not conflict_df.empty:
        print("\nCell lines with conflicting collapsed B/L labels:")
        print(
            conflict_df[
                ["cell_line", "dsmz_celldive_label", "source_labels", "collapsed_labels_observed"]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()

