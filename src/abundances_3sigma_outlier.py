
#!/usr/bin/env python3
"""
Refactored abundances module: provides AbundanceDetermination class
that integrates with myspec.SpectralAnalysisPipeline.

- No top-level execution.
- Parameterized paths and external binary names.
- Uses subprocess for external calls and ensures directories exist.
- Writes under outputs/* by default to match the project layout.
"""

from __future__ import annotations

import os
import re
import string
import subprocess
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd


def weighted_avg_and_std(values: np.ndarray) -> Tuple[float, float]:
    """Robust weighted mean/std around the median (as in your original)."""
    values = np.array(values, dtype=float)
    weights = (np.abs(values - np.median(values)) / (np.std(values) + 1e-13) + 0.25)
    weights_rounded = (np.round(weights / 0.5 + 1e-13) * 0.5) ** (-1)
    average = float(np.round(np.average(values, weights=weights_rounded), 3))
    std = float(np.sqrt(np.average((values - average) ** 2, weights=weights_rounded)))
    return average, std


def round_float(x):
    try:
        return round(float(x), 3)
    except Exception:
        return x


class AbundanceDetermination:
    def __init__(
        self,
        stellar_params: pd.DataFrame,
        atmos_dir: str = "outputs/atmospheres",
        ares_corr_dir: str = "outputs/ares_ew_corr",
        moog_cmd: str = "MOOGSILENT2019",
        damping_mode: int = 1,
        elements_info_path: str | Path = "inputs/elements.rdb",
    ) -> None:
        """
        Parameters
        ----------
        stellar_params : DataFrame
            Table with columns: star, teff, logg, feh, vtur, erteff, erlogg, erfeh, ervtur
        atmos_dir : str
            Directory with atmosphere models
        ares_corr_dir : str
            Directory with corrected ARES EW tables (+ err_ew subdir)
        moog_cmd : str
            Executable name for MOOG in silent mode
        damping_mode : int
            Damping flag (0/1/2) passed into MOOG par files
        """
        self.input_df = stellar_params.copy()
        self.atmos_dir = Path(atmos_dir)
        self.ares_corr_dir = Path(ares_corr_dir)
        self.moog_cmd = moog_cmd
        self.damping_mode = damping_mode
        self.elements_info_path = Path(elements_info_path)

        # Output roots
        self.out_root = Path("outputs/moog_abundances")
        (self.out_root / "moog_outab").mkdir(parents=True, exist_ok=True)
        (self.out_root / "moog_std_out").mkdir(parents=True, exist_ok=True)
        (self.out_root / "moog_outab_formatted").mkdir(parents=True, exist_ok=True)
        (self.out_root / "moog_abstar").mkdir(parents=True, exist_ok=True)

    # ---------- small utils ----------
    def _run(self, cmd: str):
        """Run a shell command safely; raise on failure."""
        res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if res.returncode != 0:
            raise RuntimeError(
                f"Command failed: {cmd}\nSTDOUT:\n{res.stdout.decode()}\nSTDERR:\n{res.stderr.decode()}"
            )
        return res

    def _write_text(self, path: Path, text: str):
        path.write_text(text)

    # ---------- par file writers ----------
    def abfind_par(self, damping: int):
        fout = "abfind\n"
        fout += "terminal    X11\n"
        fout += "atmosphere  1\n"
        fout += "molecules   1\n"
        fout += "lines       1\n"
        fout += "flux/int    0\n"
        fout += "plot        0\n"
        fout += f"damping     {damping}\n"
        fout += "units       0\n"
        fout += "standard_out  'standard_out'\n"
        fout += "summary_out   'outab.dat'\n"
        fout += "model_in      'out.atm'\n"
        fout += "lines_in      'lines.dat'\n"
        self._write_text(Path("abfind.par"), fout)

    def blends_par(self, damping: int):
        # Keep available if you want to extend with HFS/isotopes
        fout = "blends\n"
        fout += "terminal    X11\n"
        fout += "atmosphere  1\n"
        fout += "molecules   2\n"
        fout += "lines       1\n"
        fout += "freeform    0\n"
        fout += "flux/int    0\n"
        fout += "plot        0\n"
        fout += f"damping     {damping}\n"
        fout += "units       0\n"
        fout += "standard_out  'standard_out'\n"
        fout += "summary_out   'outab.dat'\n"
        fout += "model_in      'out.atm'\n"
        fout += "lines_in      'lines.dat'\n"
        self._write_text(Path("blends.par"), fout)

    # ---------- core routines ----------
    def abfind(self, starname: str, param: str):
        """Run MOOG abfind for star and parameter tag (orig, erteff, ...)."""
        print("MOOG abfind:", starname, param)
        self.abfind_par(self.damping_mode)

        # Prepare inputs
        self._run(f"cp {self.atmos_dir}/{starname}_{param}.atm out.atm")
        # choose EW file (normal vs err path is handled by caller)
        self._run(f"cp {self.ares_corr_dir}/{starname}_aresout.dat lines.dat")

        # Run MOOG (silent)
        self._run(f"echo abfind.par | {self.moog_cmd}")

        # Format MOOG output
        ALPHA = string.ascii_letters
        with open('outab.dat', 'r') as fin, open('outab_formated.tmp', 'w') as fout:
            for line in fin:
                if line.lstrip():
                    line = line.lstrip()
                    if not line.startswith(tuple(ALPHA)):
                        new_line = re.sub(' +','\t', line)
                        fout.write(new_line)

        abund_star = pd.read_table(
            "outab_formated.tmp",
            skiprows=1,
            usecols=(0,1,4,6,7),
            names=['wave', 'num', 'EWin', 'abund', 'delavg'],
            delimiter=r'\s+'
        )
        abund_star.to_csv(self.out_root / "moog_outab_formatted" / f"{starname}_{param}_outab_formated.dat",
                          sep='\t', index=False)

        element_list_tmp = abund_star.groupby('num')['num'].count()
        element_list = pd.DataFrame([element_list_tmp.index, element_list_tmp.values],
                                    index=['num','nlines']).transpose()

        element_info = pd.read_table(self.elements_info_path)
        element_abund = pd.merge(left=element_list, right=element_info, on='num')
        element_abund['abund'] = 0.0
        element_abund['abund_std_err'] = 0.0

        for i, num in enumerate(element_abund['num']):
            line_abundances = abund_star[abund_star['num'] == num]
            nlines = int(element_abund['nlines'][i])
            if nlines == 1:
                el_weighted_mean = line_abundances['abund']
                el_weighted_std = np.nan
            elif nlines == 2:
                el_weighted_mean = np.mean(line_abundances['abund'])
                el_weighted_std = (np.max(line_abundances['abund']) - np.min(line_abundances['abund']))/2
            else:
                removing_outliers = np.abs(line_abundances['abund'] - np.median(line_abundances['abund']))/(np.std(line_abundances['abund'])+1E-13) < 3
                line_best_abundances = line_abundances[removing_outliers]
                el_weighted_mean, el_weighted_std = weighted_avg_and_std(np.array(line_best_abundances['abund']))
                nlines = len(line_best_abundances['abund'])

            element_abund.at[i, 'nlines'] = float(nlines)
            element_abund.at[i, 'abund'] = round_float(el_weighted_mean)
            element_abund.at[i, 'abund_std_err'] = round_float(el_weighted_std)

        element_abund.to_csv(self.out_root / "moog_abstar" / f"{starname}_{param}_abund.dat",
                             sep='\t', index=False)

        # Move / clean
        os.replace("outab.dat", self.out_root / "moog_outab" / f"{starname}_{param}_outab.dat")
        os.replace("standard_out", self.out_root / "moog_std_out" / f"{starname}_{param}_std_out.dat")
        for fname in ["out.atm", "lines.dat", "abfind.par", "blends.par", "outab_formated.tmp"]:
            p = Path(fname)
            if p.exists():
                p.unlink()

        for p in Path(".").glob("fort*"):
            p.unlink()

    def abfind_elem_lt4lines(self, starname: str, param: str):
        """abfind for EW error propagation on elements with <4 lines."""
        self.abfind_par(self.damping_mode)

        self._run(f"cp {self.atmos_dir}/{starname}_{param}.atm out.atm")
        self._run(f"cp {self.ares_corr_dir}/err_ew/{starname}_aresout_ew_err.dat lines.dat")
        self._run(f"echo abfind.par | {self.moog_cmd}")

        ALPHA = string.ascii_letters
        with open('outab.dat', 'r') as fin, open('outab_formated.tmp', 'w') as fout:
            for line in fin:
                if line.lstrip():
                    line = line.lstrip()
                    if not line.startswith(tuple(ALPHA)):
                        new_line = re.sub(' +','\t', line)
                        fout.write(new_line)

        abund_star = pd.read_table(
            "outab_formated.tmp",
            skiprows=1,
            usecols=(0,1,4,6,7),
            names=['wave','num','EWin','abund','delavg'],
            delimiter=r'\s+'
        )
        abund_star.to_csv(self.out_root / "moog_outab_formatted" / f"{starname}_{param}_outab_formated_ew_err.dat",
                          sep='\t', index=False)

        element_list_tmp = abund_star.groupby('num')['num'].count()
        element_list = pd.DataFrame([element_list_tmp.index, element_list_tmp.values],
                                    index=['num','nlines']).transpose()

        element_info = pd.read_table(self.elements_info_path)
        element_abund = pd.merge(left=element_list, right=element_info, on='num')
        element_abund['abund'] = 0.0
        element_abund['abund_std_err_ew'] = 0.0

        for i, num in enumerate(element_abund['num']):
            line_abundances = abund_star[abund_star['num'] == num]
            nlines = int(element_abund['nlines'][i])
            if nlines == 1:
                el_weighted_mean = line_abundances['abund']
                el_weighted_std = np.nan
            elif nlines == 2:
                el_weighted_mean = np.mean(line_abundances['abund'])
                el_weighted_std = (np.max(line_abundances['abund']) - np.min(line_abundances['abund']))/2
            else:
                removing_outliers = np.abs(line_abundances['abund'] - np.median(line_abundances['abund']))/(np.std(line_abundances['abund'])+1E-13) < 3
                line_best_abundances = line_abundances[removing_outliers]
                el_weighted_mean, el_weighted_std = weighted_avg_and_std(np.array(line_best_abundances['abund']))
                nlines = len(line_best_abundances['abund'])

            element_abund.at[i, 'nlines'] = float(nlines)
            element_abund.at[i, 'abund'] = round_float(el_weighted_mean)
            element_abund.at[i, 'abund_std_err_ew'] = round_float(el_weighted_std)

        element_abund_selected = element_abund[['num','abund_std_err_ew']]
        element_abund_selected.to_csv(self.out_root / "moog_abstar" / f"{starname}_{param}_abund_ew_err.dat",
                                      sep='\t', index=False)

        os.replace("outab.dat", self.out_root / "moog_outab" / f"{starname}_{param}_outab_ew_err.dat")
        os.replace("standard_out", self.out_root / "moog_std_out" / f"{starname}_{param}_std_out_ew_err.dat")
        for fname in ["out.atm", "lines.dat", "abfind.par", "blends.par", "outab_formated.tmp"]:
            p = Path(fname)
            if p.exists():
                p.unlink()
        for p in Path('.').glob('fort*'):
            p.unlink()

    def abfind_err_param(self, starname: str):
        """Derive abundance sensitivity to parameter variations."""
        element_abund_error = pd.read_csv(
            self.out_root / "moog_abstar" / f"{starname}_orig_abund.dat",
            usecols=(2,8),
            skiprows=1,
            names=['element','abund_orig'],
            delimiter='\t'
        )

        element_abund_error['abund_diff_erlogg'] = 0.0
        element_abund_error['abund_diff_erteff'] = 0.0
        element_abund_error['abund_diff_erfeh'] = 0.0
        element_abund_error['abund_diff_ervtur'] = 0.0

        parameters = ['erlogg', 'erteff', 'erfeh', 'ervtur']
        for parameter in parameters:
            self.abfind(starname, parameter)
            slected_columns = pd.read_csv(
                self.out_root / "moog_abstar" / f"{starname}_{parameter}_abund.dat",
                usecols=(2,8),
                skiprows=1,
                names=['element', f'abund_{parameter}'],
                delimiter='\t'
            )
            merged = element_abund_error.merge(slected_columns, on='element')
            merged[f'abund_diff_{parameter}'] = np.abs(merged['abund_orig'] - merged[f'abund_{parameter}'])
            element_abund_error[f'abund_diff_{parameter}'] = merged[f'abund_diff_{parameter}']

        element_abund_error['abund_param_err'] = (
            element_abund_error['abund_diff_erlogg']**2 +
            element_abund_error['abund_diff_erteff']**2 +
            element_abund_error['abund_diff_erfeh']**2 +
            element_abund_error['abund_diff_ervtur']**2
        )**0.5
        element_abund_param_error = element_abund_error[['element','abund_param_err']]

        element_abund = pd.read_table(self.out_root / "moog_abstar" / f"{starname}_orig_abund.dat", delimiter='\t')
        element_abund_all = element_abund.merge(element_abund_param_error, on='element')
        element_abund_all.to_csv(self.out_root / "moog_abstar" / f"{starname}_abund_param_err.dat",
                                 sep='\t', index=False)

        # clean parameter intermediates (mirror original behavior)
        for suffix in ['erlogg','erteff','erfeh','ervtur']:
            for sub in ['moog_abstar','moog_outab_formatted','moog_outab','moog_std_out']:
                for p in (self.out_root / sub).glob(f'*{suffix}*'):
                    p.unlink(missing_ok=True)

    def abund_err(self, starname: str):
        self.abfind_elem_lt4lines(starname, 'orig')

        element_abund = pd.read_table(self.out_root / "moog_abstar" / f"{starname}_abund_param_err.dat", delimiter='\t')
        element_abund_ew_err = pd.read_table(self.out_root / "moog_abstar" / f"{starname}_orig_abund_ew_err.dat", delimiter='\t')
        element_abund_all = element_abund.merge(element_abund_ew_err, on='num', how='outer')
        element_abund_all['total_err'] = 0

        mask_has_ew = element_abund_all['abund_std_err_ew'] > 0
        element_abund_all.loc[mask_has_ew, 'total_err'] = (element_abund_all['abund_std_err_ew']**2 + element_abund_all['abund_param_err']**2)**0.5
        mask_nan = element_abund_all['abund_std_err_ew'] != element_abund_all['abund_std_err_ew']
        element_abund_all.loc[mask_nan, 'total_err'] = (element_abund_all['abund_std_err']**2 + element_abund_all['abund_param_err']**2)**0.5

        element_abund_all['abund_param_err'] = np.round(element_abund_all['abund_param_err'], 3)
        element_abund_all['total_err'] = np.round(element_abund_all['total_err'], 3)
        element_abund_all.to_csv(self.out_root / "moog_abstar" / f"{starname}_abund_err.dat", sep='\t', index=False)

    # ---------- public entry ----------
    def calculate_abundances(self):
        """
        Main entry: derive abundances + errors for all stars in self.input_df.
        Requires the `outputs/ares_ew_corr` and `outputs/atmospheres` steps to be completed.
        """
        for i in range(len(self.input_df)):
            star = str(self.input_df.star[i])
            print("Processing:", star)
            # main abundances
            self.abfind(star, 'orig')
            # parameter sensitivities
            self.abfind_err_param(star)
            # EW-based error and combine
            self.abund_err(star)
            # remove intermediate error files as in original
            for pat in [f"{star}_abund_param_err.dat", f"{star}_orig_abund.dat", f"{star}_orig_abund_ew_err.dat"]:
                for sub in ["moog_abstar"]:
                    p = self.out_root / sub / pat
                    if p.exists():
                        p.unlink(missing_ok=True)
