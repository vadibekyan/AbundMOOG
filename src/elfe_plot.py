
#!/usr/bin/python
"""
ELFE plotting utilities (refactored):
- No top-level execution.
- Parameterized inputs (files, target_name).
- Safe to import and call from myspec.py pipeline.
- Produces PDFs in a specified output directory.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def N_columns_rows_for_plot(N_elements: int) -> Tuple[int, int]:
    n_rows = 0
    n_columns = 0
    if N_elements == 4:
        n_rows, n_columns = 2, 2
    elif 5 <= N_elements <= 6:
        n_rows, n_columns = 3, 2
    elif 7 <= N_elements <= 9:
        n_rows, n_columns = 3, 3
    elif 10 <= N_elements <= 12:
        n_rows, n_columns = 4, 3
    elif 13 <= N_elements <= 15:
        n_rows, n_columns = 5, 3
    elif N_elements == 16:
        n_rows, n_columns = 4, 4
    elif 17 <= N_elements <= 20:
        n_rows, n_columns = 5, 4
    elif 21 <= N_elements <= 24:
        n_rows, n_columns = 6, 4
    elif 22 <= N_elements <= 28:
        n_rows, n_columns = 7, 4
    return n_rows, n_columns

def tick_parameters(n_rows: int, n_columns: int):
    top_left = ()
    for i in range(n_rows-1):
        top_left += +i*n_columns,
    top_left += 100,

    top_right = ()
    for j in range(n_columns-1):
        for i in range(n_rows-1):
            top_right += (1+i*n_columns)+(j),
    top_right += 100,

    bottom_left = ()
    bottom_left += (n_rows-1)*n_columns,
    bottom_left += 100,

    bottom_right = ()
    for j in range(n_columns-1):
        bottom_right += (n_rows-1)*n_columns + j+1,
    bottom_right += 100,
    return bottom_left, top_left, top_right, bottom_right

class ElfePlotter:
    def __init__(
        self,
        abund_star_table: Path | str = "all_abundances_ordered_star.rdb",
        elements_to_plot: Path | str = "elements_to_plot.rdb",
        harps_catalog: Path | str = "abundances_harps1111_o_c.rdb",
        target_name: str = "",
        outdir: Path | str = "plots",
    ) -> None:
        self.abund_star_table = Path(abund_star_table)
        self.elements_to_plot = Path(elements_to_plot)
        self.harps_catalog = Path(harps_catalog)
        self.target_name = target_name
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)

        # Load once
        self.abundances = pd.read_table(self.abund_star_table)
        self.elements_df = pd.read_table(self.elements_to_plot)
        self.harps_sample = pd.read_table(self.harps_catalog, sep='\t')
        self.harps_solar_analogs = self.harps_sample[
            (self.harps_sample['teff'] > 5277) & (self.harps_sample['teff'] < 6277) & (self.harps_sample['feh'] < 1.4)
        ]

    def _apply_style(self):
        sns.set_style("white")
        sns.set_style("ticks")
        sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    def elfe_feh(self) -> Path:
        self._apply_style()
        n_rows, n_columns = N_columns_rows_for_plot(len(self.elements_df.element))

        # Collect valid [X/Fe] values across all elements to set global y-limits
        all_vals = []
        for element in self.elements_df.element:
            if element not in self.abundances.columns:
                continue
            abundances_df = self.abundances[pd.notna(self.abundances.get(f'{element}_rel_sun'))]
            if abundances_df.empty:
                continue
            elfe = abundances_df[f'{element}_rel_sun'] - abundances_df['feh'] + 0.02
            elfe = elfe.replace([np.inf, -np.inf], np.nan).dropna()
            if not elfe.empty:
                all_vals.extend(elfe.tolist())

        # Fallback limits if we didn't get any data
        if not all_vals:
            min_all_elfe, max_all_elfe = -0.5, 0.8
        else:
            min_all_elfe = np.floor(np.nanmin(all_vals) * 10) / 10 - 0.1
            max_all_elfe = np.ceil(np.nanmax(all_vals) * 10) / 10

        plt.figure(1)
        fig = plt.gcf()
        fig.set_size_inches(3*n_columns, 3*n_rows)
        fig.subplots_adjust(left=0.13, right=0.98, bottom=0.09, top=0.98, wspace=0.0, hspace=0.0)

        for i, element in enumerate(self.elements_df.element):
            # Safely get comparison sample column if it exists
            harps_col_ok = element in self.harps_solar_analogs.columns
            if harps_col_ok:
                harps_best = self.harps_solar_analogs[
                    (self.harps_solar_analogs[element] < 1) &
                    (self.harps_solar_analogs[element] > -2)
                ]
            else:
                harps_best = self.harps_solar_analogs.iloc[0:0]  # empty

            # Abundance rows for this element
            if element not in self.abundances.columns or f'{element}_rel_sun' not in self.abundances.columns:
                abundances_df = self.abundances.iloc[0:0]  # empty
            else:
                abundances_df = self.abundances[pd.notna(self.abundances[f'{element}_rel_sun'])]

            # Compute ELFE (+ errors) safely
            if not abundances_df.empty:
                elfe = abundances_df[f'{element}_rel_sun'] - abundances_df['feh'] + 0.02
                elfe_err = ((abundances_df.get(f'{element}_err', 0)**2) + (abundances_df.get('erfeh', 0)**2))**0.5
                giants = abundances_df[abundances_df.teff < 4300]
                dwarfs = abundances_df[abundances_df.logg >= 0.5]
                elfe_giants = giants[f'{element}_rel_sun'] - giants['feh'] + 0.02
                elfe_err_giants = ((giants.get(f'{element}_err', 0)**2) + (giants.get('erfeh', 0)**2))**0.5
                elfe_dwarfs = dwarfs[f'{element}_rel_sun'] - dwarfs['feh'] + 0.02
                elfe_err_dwarfs = ((dwarfs.get(f'{element}_err', 0)**2) + (dwarfs.get('erfeh', 0)**2))**0.5
                target_abund = abundances_df[abundances_df.star == self.target_name]
                if not target_abund.empty:
                    elfe_target_abund = target_abund[f'{element}_rel_sun'] - target_abund['feh'] + 0.02
                    elfe_err_target_abund = ((target_abund.get(f'{element}_err', 0)**2) + (target_abund.get('erfeh', 0)**2))**0.5
            else:
                elfe = elfe_err = giants = dwarfs = elfe_giants = elfe_err_giants = elfe_dwarfs = elfe_err_dwarfs = target_abund = None

            ax = plt.subplot(n_rows, n_columns, i+1)
            ax.text(.9, .85, f'{element}', ha='right', transform=ax.transAxes, fontsize=14)

            # Plot comparison sample if available
            if harps_col_ok and not harps_best.empty:
                plt.scatter(harps_best.feh, (harps_best[element]-harps_best.feh), s=20, alpha=0.5, edgecolors='none')

            # Plot dwarfs/giants if we have data
            if dwarfs is not None and not dwarfs.empty:
                plt.scatter(dwarfs['feh'], elfe_dwarfs, linewidth=1, alpha=0.6, s=70)
                plt.errorbar(dwarfs['feh'], elfe_dwarfs, yerr=elfe_err_dwarfs, xerr=dwarfs.get('erfeh', 0), fmt='o', markersize=1, linewidth=0.5, alpha=0.6)
            if giants is not None and not giants.empty:
                plt.scatter(giants['feh'], elfe_giants, linewidth=1, alpha=0.6, s=90)
                plt.errorbar(giants['feh'], elfe_giants, yerr=elfe_err_giants, xerr=giants.get('erfeh', 0), fmt='o', markersize=1, linewidth=0.5, alpha=0.6)

            # Highlight target if present
            if target_abund is not None and not target_abund.empty:
                plt.scatter(target_abund['feh'], elfe_target_abund, linewidth=1, alpha=0.9, s=150)
                plt.errorbar(target_abund['feh'], elfe_target_abund, yerr=elfe_err_target_abund, xerr=target_abund.get('erfeh', 0), linewidth=2, alpha=0.99, ms=15)

            plt.axhline(y=0, linestyle='--', linewidth=1.)
            plt.axvline(x=0, linestyle='--', linewidth=1.)
            #plt.ylim(min_all_elfe, max_all_elfe)
            plt.ylim(-0.3, 0.6)
            plt.xlim(-0.9, 0.5)
            plt.tick_params(labelsize=13)

            bottom_left, top_left, top_right, bottom_right = tick_parameters(n_rows, n_columns)
            if i in top_right:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=False)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=False, labelsize=13)
            if i in top_left:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=True, labelsize=13)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=False, labelsize=14)
            if i in bottom_right:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=False, labelsize=14)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=True, labelsize=12.5)
            if i in bottom_left:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=True, labelsize=13)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=True, labelsize=12.5)

        fig.text(0.55, 0.01, '[Fe/H]', ha='center', fontsize=17)
        fig.text(0.01, 0.5, '[X/Fe]', va='center', rotation='vertical', fontsize=17)
        out = self.outdir / 'elfe_feh.pdf'
        plt.savefig(out, format='pdf')
        plt.clf()
        return out


    def elfe_teff(self) -> Path:
        self._apply_style()
        n_rows, n_columns = N_columns_rows_for_plot(len(self.elements_df.element))

        plt.figure(1)
        fig = plt.gcf()
        fig.set_size_inches(3*n_columns, 3*n_rows)
        fig.subplots_adjust(left=0.1, right=0.98, bottom = 0.07, top = 0.98, wspace=0.0, hspace=0.0)

        for i, element in enumerate(self.elements_df.element):
            harps_best = self.harps_sample[(self.harps_sample[element] < 1) & (self.harps_sample[element] > -2)]
            abundances_df = self.abundances[(self.abundances[f'{element}_rel_sun'] > -2)]
            elfe = abundances_df[f'{element}_rel_sun'] - abundances_df['feh'] + 0.02
            elfe_err = (abundances_df[f'{element}_err']**2 + abundances_df['erfeh']**2)**0.5

            giants = abundances_df[abundances_df.logg < 3.5]
            elfe_giants = giants[f'{element}_rel_sun'] - giants['feh'] + 0.02
            elfe_err_giants = (giants[f'{element}_err']**2 + giants['erfeh']**2)**0.5

            dwarfs = abundances_df[abundances_df.logg >= 3.5]
            elfe_dwarfs = dwarfs[f'{element}_rel_sun'] - dwarfs['feh'] + 0.02
            elfe_err_dwarfs = (dwarfs[f'{element}_err']**2 + dwarfs['erfeh']**2)**0.5

            ax = plt.subplot(n_rows, n_columns, i+1)
            ax.text(.9, .85, f'{element}', ha='right', transform=ax.transAxes, fontsize=14)
            plt.scatter(harps_best.teff, (harps_best[element]-harps_best.feh), c='grey', s=20, alpha=0.5, edgecolors='none')
            plt.scatter(dwarfs['teff'], elfe_dwarfs, linewidth=1, color='black', alpha=0.6, s=70)
            plt.errorbar(dwarfs['teff'], elfe_dwarfs, yerr=elfe_err_dwarfs, xerr=dwarfs['erfeh'], color='black', fmt='o', markersize=1, linewidth=0.5, alpha=0.6)
            plt.scatter(giants['teff'], elfe_giants, linewidth=1, color='red', alpha=0.6, s=90)
            plt.errorbar(giants['teff'], elfe_giants, yerr=elfe_err_giants, xerr=giants['erfeh'], color='red', fmt='o', markersize=1, linewidth=0.5, alpha=0.6)
            plt.axhline(y=0, color='blue', linestyle='--', linewidth=1.)
            plt.axvline(x=5777, color='blue', linestyle='--', linewidth=1.)
            plt.ylim(-0.5, 0.8)
            plt.xlim(4400, 6800)
            plt.xticks(np.arange(4500, 6600, 500))
            plt.tick_params(labelsize=13)

            bottom_left, top_left, top_right, bottom_right = tick_parameters(n_rows, n_columns)
            if i in top_right:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=False)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=False, labelsize=16)
            if i in top_left:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=True, labelsize=13)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=False, labelsize=16)
            if i in bottom_right:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=False, labelsize=16)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=True, labelsize=12.5)
            if i in bottom_left:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=True, labelsize=13)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=True, labelsize=12.5)

        fig.text(0.55, 0.01, '$T_{eff}$', ha='center', fontsize=17)
        fig.text(0.02, 0.5, '[X/Fe]', va='center', rotation='vertical', fontsize=17)
        out = self.outdir / 'elfe_teff.pdf'
        plt.savefig(out, format='pdf')
        plt.clf()
        return out

    def elfe_logg(self) -> Path:
        self._apply_style()
        n_rows, n_columns = N_columns_rows_for_plot(len(self.elements_df.element))

        plt.figure(1)
        fig = plt.gcf()
        fig.set_size_inches(3*n_columns, 3*n_rows)
        fig.subplots_adjust(left=0.1, right=0.98, bottom = 0.07, top = 0.98, wspace=0.0, hspace=0.0)

        for i, element in enumerate(self.elements_df.element):
            harps_best = self.harps_sample[(self.harps_sample[element] < 1) & (self.harps_sample[element] > -2)]
            abundances_df = self.abundances[(self.abundances[f'{element}_rel_sun'] > -2)]
            elfe = abundances_df[f'{element}_rel_sun'] - abundances_df['feh'] + 0.02
            elfe_err = (abundances_df[f'{element}_err']**2 + abundances_df['erfeh']**2)**0.5

            giants = abundances_df[abundances_df.logg < 3.5]
            elfe_giants = giants[f'{element}_rel_sun'] - giants['feh'] + 0.02
            elfe_err_giants = (giants[f'{element}_err']**2 + giants['erfeh']**2)**0.5

            dwarfs = abundances_df[abundances_df.logg >= 3.5]
            elfe_dwarfs = dwarfs[f'{element}_rel_sun'] - dwarfs['feh'] + 0.02
            elfe_err_dwarfs = (dwarfs[f'{element}_err']**2 + dwarfs['erfeh']**2)**0.5

            plt.subplot(n_rows, n_columns, i+1)
            plt.text(.8, .85, f'{element}', ha='right', transform=plt.gca().transAxes, fontsize=14)
            plt.scatter(harps_best.logg, (harps_best[element]-harps_best.feh), c='grey', s=20, alpha=0.5, edgecolors='none')
            plt.scatter(dwarfs['logg'], elfe_dwarfs, linewidth=1, color='black', alpha=0.6, s=70)
            plt.errorbar(dwarfs['logg'], elfe_dwarfs, yerr=elfe_err_dwarfs, xerr=dwarfs['erfeh'], color='black', fmt='o', markersize=1, linewidth=0.5, alpha=0.6)
            plt.scatter(giants['logg'], elfe_giants, linewidth=1, color='red', alpha=0.6, s=90)
            plt.errorbar(giants['logg'], elfe_giants, yerr=elfe_err_giants, xerr=giants['erfeh'], color='red', fmt='o', markersize=1, linewidth=0.5, alpha=0.6)
            plt.axhline(y=0, color='blue', linestyle='--', linewidth=1.)
            plt.axvline(x=4.44, color='blue', linestyle='--', linewidth=1.)
            plt.ylim(-0.5, 0.8)
            plt.xlim(3.25, 4.75)
            plt.tick_params(labelsize=13)

            bottom_left, top_left, top_right, bottom_right = tick_parameters(n_rows, n_columns)
            if i in top_right:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=False)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=False)
            if i in top_left:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=True, labelsize=13)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=False)
            if i in bottom_right:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=False)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=True, labelsize=12)
            if i in bottom_left:
                plt.tick_params(axis='y', which='both', right=True, left=True, labelleft=True, labelsize=13)
                plt.tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=True, labelsize=12.0)

        out = self.outdir / 'elfe_logg.pdf'
        plt.savefig(out, format='pdf')
        plt.gcf().text(0.55, 0.01, '$\\log g$', ha='center', fontsize=17)
        plt.gcf().text(0.02, 0.5, '[X/Fe]', va='center', rotation='vertical', fontsize=17)
        plt.clf()
        return out

    def hr(self) -> Path:
        self._apply_style()
        plt.figure(1)
        plt.scatter(self.abundances['teff'], self.abundances['logg'], linewidth=1, color='black')
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        plt.ylabel('log g')
        plt.xlabel('Teff')
        out = self.outdir / 'HR.pdf'
        plt.savefig(out, format='pdf')
        plt.clf()
        return out

    def plot_all(self) -> dict:
        return {
            'elfe_feh': self.elfe_feh(),
            'elfe_teff': self.elfe_teff(),
            'elfe_logg': self.elfe_logg(),
            'HR': self.hr(),
        }
