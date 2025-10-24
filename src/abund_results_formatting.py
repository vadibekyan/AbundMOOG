#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Tuple
import numpy as np
import pandas as pd

class AbundanceResultsFormatter:
    """
    Build two formatted tables from MOOG outputs:
      - outputs/analysis/all_abundances_ordered_element.rdb  (elements in rows; per-star columns)
      - outputs/analysis/all_abundances_ordered_star.rdb     (stars in rows; per-element columns)

    Defaults assume the project layout:
      inputs/input_param_error.rdb
      inputs/elements_formatting.rdb
      outputs/moog_abundances/moog_abstar/<star>_abund_err.dat
    """

    def __init__(
        self,
        input_starlist: str | Path = "inputs/input_param_error.rdb",
        element_file:   str | Path = "inputs/elements_formatting.rdb",
        moog_root:      str | Path = "outputs/moog_abundances/moog_abstar",
    ) -> None:
        self.input_starlist = Path(input_starlist)
        self.element_file   = Path(element_file)
        self.moog_root      = Path(moog_root)

    def _read_starlist(self) -> pd.DataFrame:
        # starlist: space-separated; second header line skipped (as in your original)
        return pd.read_csv(self.input_starlist, sep=r"\s+", skiprows=[1])

    def _read_elements(self) -> pd.DataFrame:
        # elements_formatting.rdb should include at least ['element'] and optionally ['sun_vesta']
        df = pd.read_table(self.element_file)
        if "element" not in df.columns:
            raise ValueError("elements_formatting.rdb must contain a column named 'element'.")
        if "sun_vesta" not in df.columns:
            df["sun_vesta"] = 0.0
        return df

    def format_columns(self, output_file: str | Path = "outputs/analysis/all_abundances_ordered_element.rdb") -> Path:
        input_data   = self._read_starlist()
        element_info = self._read_elements().copy()

        for star in input_data["star"]:
            star_abund_file = self.moog_root / f"{star}_abund_err.dat"
            if not star_abund_file.exists():
                # still create columns so every star has the same schema
                element_info[f"{star}_abund"]         = np.nan
                element_info[f"{star}_abund_err"]     = np.nan
                element_info[f"{star}_abund_rel_sun"] = np.nan
                continue

            star_abund = pd.read_table(star_abund_file)[["element", "abund", "total_err"]]
            merged = element_info[["element", "sun_vesta"]].merge(star_abund, on="element", how="left")
            element_info[f"{star}_abund"]         = merged["abund"]
            element_info[f"{star}_abund_err"]     = merged["total_err"]
            element_info[f"{star}_abund_rel_sun"] = merged["abund"] - merged["sun_vesta"]

        out = Path(output_file)
        out.parent.mkdir(parents=True, exist_ok=True)
        element_info.to_csv(out, sep="\t", index=False, float_format="%.3f")
        return out

    def format_rows(
        self,
        column_file: str | Path = "outputs/analysis/all_abundances_ordered_element.rdb",
        output_file: str | Path = "outputs/analysis/all_abundances_ordered_star.rdb",
    ) -> Path:
        input_data = self._read_starlist()
        abund_col  = pd.read_table(column_file)

        # Build per-element columns for every star (lookup by element row)
        for element in abund_col["element"]:
            # ensure columns exist
            if element not in input_data.columns:
                input_data[f"{element}"]         = np.nan
                input_data[f"{element}_err"]     = np.nan
                input_data[f"{element}_rel_sun"] = np.nan

            # row index for this element in abund_col
            idx = abund_col.index[abund_col["element"] == element][0]
            for i, star in enumerate(input_data["star"]):
                input_data.at[i, f"{element}"]         = abund_col.get(f"{star}_abund", pd.Series([np.nan]*len(abund_col))).iat[idx]
                input_data.at[i, f"{element}_err"]     = abund_col.get(f"{star}_abund_err", pd.Series([np.nan]*len(abund_col))).iat[idx]
                input_data.at[i, f"{element}_rel_sun"] = abund_col.get(f"{star}_abund_rel_sun", pd.Series([np.nan]*len(abund_col))).iat[idx]

        out = Path(output_file)
        out.parent.mkdir(parents=True, exist_ok=True)
        input_data.to_csv(out, sep="\t", index=False, float_format="%.3f")
        return out

    def run(
        self,
        output_column_file: str | Path = "outputs/analysis/all_abundances_ordered_element.rdb",
        output_line_file:   str | Path = "outputs/analysis/all_abundances_ordered_star.rdb",
    ) -> tuple[Path, Path]:
        col = self.format_columns(output_column_file)
        row = self.format_rows(col, output_line_file)
        return col, row

if __name__ == "__main__":
    # Optional: quick self-test (won't run on import)
    afr = AbundanceResultsFormatter()
    afr.run()
