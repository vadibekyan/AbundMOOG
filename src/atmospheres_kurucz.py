#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refactored AtmosphereModeling:
- Writes all .atm files into outputs/atmospheres
- No top-level execution
- No shell mv; uses pathlib
"""

import os
from pathlib import Path
import pandas as pd


class AtmosphereModeling:
    def __init__(self, input_params: pd.DataFrame, atmos_out_dir: str | Path = "outputs/atmospheres"):
        self.input_df = input_params.copy()
        self.atmos_out_dir = Path(atmos_out_dir)
        self.atmos_out_dir.mkdir(parents=True, exist_ok=True)

    def _run_intermod_transform(self, teff: float, logg: float, feh: float, vtur: float, dest_name: str):
        """
        Run external tools to produce out.atm, then move it to outputs/atmospheres/<dest_name>.atm
        """
        # External binaries must be available in PATH:
        # intermod_kurucz reads Teff, logg, [M/H] from stdin, writes mod* files
        os.system(f'echo {teff} {logg} {feh} | intermod_kurucz')
        # transform_kurucz reads vturb from stdin, writes out.atm
        os.system(f'echo {vtur} | transform_kurucz')

        dest = self.atmos_out_dir / f"{dest_name}.atm"
        Path("out.atm").replace(dest)

    def atmo_kurucz(self, starname, teff, logg, feh, vtur, erteff, erlogg, erfeh, ervtur):
        """Create the five atmosphere models directly into outputs/atmospheres/."""
        # orig
        self._run_intermod_transform(teff, logg, feh, vtur, f"{starname}_orig")
        # erfeh
        self._run_intermod_transform(teff, logg, feh + erfeh, vtur, f"{starname}_erfeh")
        # erteff
        self._run_intermod_transform(teff + erteff, logg, feh, vtur, f"{starname}_erteff")
        # erlogg
        self._run_intermod_transform(teff, logg + erlogg, feh, vtur, f"{starname}_erlogg")
        # ervtur
        self._run_intermod_transform(teff, logg, feh, vtur + ervtur, f"{starname}_ervtur")

        # clean intermod artifacts if any
        os.system('rm -f mod*')

    def model_atmospheres(self):
        """Main loop."""
        starlist = self.input_df
        for i in range(len(starlist)):
            self.atmo_kurucz(
                starname=str(starlist.star[i]),
                teff=float(starlist.teff[i]),
                logg=float(starlist.logg[i]),
                feh=float(starlist.feh[i]),
                vtur=float(starlist.vtur[i]),
                erteff=float(starlist.erteff[i]),
                erlogg=float(starlist.erlogg[i]),
                erfeh=float(starlist.erfeh[i]),
                ervtur=float(starlist.ervtur[i]),
            )
