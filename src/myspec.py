#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
myspec.py — pipeline entrypoint (config-driven, import-safe)

Run from project root:
    python -m src.myspec --config config.json --target-name "Target-123"
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from .ares_ew import EWMeasurements
from .atmospheres_kurucz import AtmosphereModeling
from .abundances_3sigma_outlier import AbundanceDetermination
from .abund_results_formatting import AbundanceResultsFormatter
from .elfe_plot import ElfePlotter


class SpectralAnalysisPipeline:
    def __init__(self, stellar_params: pd.DataFrame):
        self.stellar_params = stellar_params
        self.ew_measurements = EWMeasurements(stellar_params)
        self.atmosphere_modeling = None
        self.abundance_determination = None

    def measure_EWs(self):
        self.ew_measurements.ew_measurements()
        print("Equivalent Widths (EWs) measured.")

    def model_atmospheres(self):
        self.atmosphere_modeling = AtmosphereModeling(self.stellar_params)
        self.atmosphere_modeling.model_atmospheres()
        print("Atmosphere models created.")

    def determine_abundances(self):
        self.abundance_determination = AbundanceDetermination(self.stellar_params)
        self.abundance_determination.calculate_abundances()
        print("Abundances determined.")

    def format_abundance_outputs(
        self,
        element_file: str = "inputs/elements_formatting.rdb",
        out_column: str = "outputs/analysis/all_abundances_ordered_element.rdb",
        out_line: str = "outputs/analysis/all_abundances_ordered_star.rdb",
    ):
        """Format abundances into element-wise and star-wise tables."""
        formatter = AbundanceResultsFormatter(
            input_starlist="inputs/input_param_error.rdb",
            element_file=element_file,
            moog_root="outputs/moog_abundances/moog_abstar",
        )
        col_path, row_path = formatter.run(out_column, out_line)
        print(f"Abundance tables written to: {col_path} and {row_path}")

    def plot_abundances(
        self,
        elements_to_plot: str = "inputs/elements_to_plot.rdb",
        harps_catalog: str = "inputs/abundances_harps1111_o_c.rdb",
        target_name: str = "",
        outdir: str = "outputs/plots",
    ):
        plotter = ElfePlotter(
            abund_star_table="outputs/analysis/all_abundances_ordered_star.rdb",
            elements_to_plot=elements_to_plot,
            harps_catalog=harps_catalog,
            target_name=target_name,
            outdir=outdir,
        )
        outputs = plotter.plot_all()
        print("Plot files written:", outputs)

    def run_pipeline(self):
        self.measure_EWs()
        self.model_atmospheres()
        self.determine_abundances()
        self.format_abundance_outputs()
        self.plot_abundances()
        print("Spectral analysis pipeline completed.")


def _load_config(path: str | None) -> dict:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Config file not found: {p}")
    with open(p, "r") as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description="Run myspec pipeline (config-driven).")
    parser.add_argument("--config", type=str, default=None, help="Path to config.json")
    parser.add_argument("--target-name", type=str, default="", help="Name of target star to highlight in plots")
    args = parser.parse_args()

    cfg = _load_config(args.config)

    # Load stellar params from inputs (respect layout); allow override via config if provided
    starlist_path = cfg.get("formatting", {}).get("input_starlist", "inputs/input_param_error.rdb")
    stellar_params = pd.read_csv(starlist_path, sep=r"\s+", skiprows=[1])

    pipeline = SpectralAnalysisPipeline(stellar_params)

    # Run full pipeline
    pipeline.run_pipeline()

    # Optionally re-run plots with overrides
    plot_cfg = cfg.get("plot", {})
    elements_to_plot = plot_cfg.get("elements_to_plot", "inputs/elements_to_plot.rdb")
    harps_catalog = plot_cfg.get("harps_catalog", "inputs/abundances_harps1111_o_c.rdb")
    outdir = plot_cfg.get("outdir", "outputs/plots")
    target_name = args.target_name or plot_cfg.get("target_name", "")

    pipeline.plot_abundances(
        elements_to_plot=elements_to_plot,
        harps_catalog=harps_catalog,
        target_name=target_name,
        outdir=outdir,
    )


if __name__ == "__main__":
    main()
