"""Export AffiTox Pearson r bootstrap CIs for SI Table tab:si-pearson-ci.

Reads analysis/tables/correlations/correlations_with_ci_perm.csv and writes
LaTeX tabular rows plus a tidy CSV under the same directory.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "analysis" / "tables" / "correlations" / "correlations_with_ci_perm.csv"
OUT_CSV = ROOT / "analysis" / "tables" / "correlations" / "pearson_ci_for_si.csv"
OUT_TEX = ROOT / "analysis" / "tables" / "correlations" / "pearson_ci_si_rows.tex"

ORDER = ["boltz2", "dynamicbind", "gnina", "plapt", "qvina"]
LABELS = {
    "boltz2": "Boltz-2",
    "dynamicbind": "DynamicBind",
    "gnina": "GNINA 1.3",
    "plapt": "PLAPT",
    "qvina": "QVina2",
}


def _num(x: float) -> str:
    return f"$-${abs(x):.2f}" if x < 0 else f"{x:.2f}"


def _val(x: float) -> str:
    return f"$-${abs(x):.3f}" if x < 0 else f"{x:.3f}"


def main() -> None:
    df = pd.read_csv(SRC)
    rows = []
    tex_lines = []
    for tgt in sorted(df["target"].unique()):
        cells = [tgt.upper()]
        g = df[df["target"] == tgt]
        for mid in ORDER:
            r = g[g["method_id"] == mid].iloc[0]
            star = r["pearson_p_bh"] < 0.05
            lo, hi = float(r["pearson_ci_low"]), float(r["pearson_ci_high"])
            val = float(r["pearson_r"])
            cell = f"{_val(val)}$^*$ [{_num(lo)},{_num(hi)}]" if star else f"{_val(val)} [{_num(lo)},{_num(hi)}]"
            cells.append(cell)
            rows.append(
                {
                    "target": tgt.upper(),
                    "method": LABELS[mid],
                    "n_points": int(r["n_points"]),
                    "pearson_r": val,
                    "pearson_ci_low": lo,
                    "pearson_ci_high": hi,
                    "pearson_p_bh": float(r["pearson_p_bh"]),
                    "bh_sig": bool(star),
                }
            )
        tex_lines.append(" & ".join(cells) + r" \\")

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    OUT_TEX.write_text("\n".join(tex_lines) + "\n")
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_TEX}")


if __name__ == "__main__":
    main()
