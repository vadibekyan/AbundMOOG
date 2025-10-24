
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ARES EW pipeline (refactored, import-safe)

- Reads line list from inputs/linelist_damp.rdb
- Looks for spectra under inputs/spectra/<star>.fits (fallback: inputs/<star>.fits)
- Runs ARES in the current working directory (expects ARES in PATH)
- Writes results under outputs/ares_* and outputs/ares_ew_corr
- No top-level execution (safe to import)
"""

import os
import decimal
from pathlib import Path
import numpy as np
import pandas as pd


class EWMeasurements:
    def __init__(self, input_params: pd.DataFrame,
                 inputs_dir: str | Path = "inputs",
                 outputs_dir: str | Path = "outputs"):
        self.input_df = input_params.copy()
        self.inputs_dir = Path(inputs_dir)
        self.outputs_dir = Path(outputs_dir)

        # Output directories
        self.dir_log = self.outputs_dir / "ares_log"
        self.dir_out = self.outputs_dir / "ares_out"
        self.dir_ew = self.outputs_dir / "ares_ew"
        self.dir_ew_err = self.dir_ew / "err_ew"
        self.dir_corr = self.outputs_dir / "ares_ew_corr"

        for d in [self.dir_log, self.dir_out, self.dir_ew, self.dir_ew_err, self.dir_corr]:
            d.mkdir(parents=True, exist_ok=True)

        # Inputs
        self.linelist_path = self.inputs_dir / "linelist_damp.rdb"
        if not self.linelist_path.exists():
            raise FileNotFoundError(f"Line list not found: {self.linelist_path}")

    @staticmethod
    def round_up0(i):
        """Round to 2 decimals with HALF_UP to match original behavior."""
        rounded = decimal.Decimal(str(i)).quantize(decimal.Decimal('1.11'), rounding=decimal.ROUND_HALF_UP)
        return float(rounded)

    @staticmethod
    def _read_snr_from_log(logfile: Path) -> int:
        with open(logfile, 'r') as lines:
            for line in lines:
                if line.startswith('S/N'):
                    parts = line.strip('\n').split(':')
                    return int(float(parts[1]))
        raise RuntimeError(f"Could not find S/N in {logfile}")

    def _write_mine_opt(self, spectra_name: str):
        """Create mine.opt in CWD for ARES. Prefers inputs/spectra/<star>.fits; fallback inputs/<star>.fits"""
        spec_in_spectra = self.inputs_dir / "spectra" / f"{spectra_name}.fits"
        spec_in_inputs = self.inputs_dir / f"{spectra_name}.fits"
        if spec_in_spectra.exists():
            specfits_line = f"specfits='./{spec_in_spectra.as_posix()}'\n"
        elif spec_in_inputs.exists():
            specfits_line = f"specfits='./{spec_in_inputs.as_posix()}'\n"
        else:
            # Leave relative name; user may have it in CWD
            specfits_line = f"specfits='./{spectra_name}.fits'\n"

        fout  = specfits_line
        fout += "readlinedat= 'cdo.dat'\n"
        fout += "fileout='aresout.dat'\n"
        fout += "lambdai= 3000\n"
        fout += "lambdaf=8000\n"
        fout += "smoothder=14\n"
        fout += "space=3.0\n"
        fout += "rejt= 3;5764,5766,6047,6052,6068,6076\n"
        fout += "lineresol=0.1\n"
        fout += "miniline=2\n"
        fout += "plots_flag=0\n"
        fout += "rvmask='3,6021.8,6024.06,6027.06,6024.06,20'\n"
        with open('mine.opt', 'w') as f:
            f.writelines(fout)

    def _write_cdo_from_linelist(self):
        """Create cdo.dat with the first column (WL) from inputs/linelist_damp.rdb, sorted, no rounding loss."""
        linelist = pd.read_csv(
            self.linelist_path,
            skiprows=2,
            names=['WL', 'EP', 'loggf', 'element', 'num', 'damp', 'EWsun'],
            sep=r'\s+',
            engine='python',
            dtype={'WL': float},
        ).sort_values(['WL'])  # keep your original behavior
        # ARES expects a simple column of wavelengths; write ALL rows in order
        with open('cdo.dat', 'w') as f:
            for wl in linelist['WL']:
                f.write(f"{wl:.2f}\n")  # format with 2 decimals only for the file, not for internal math
        # keep a copy for downstream methods
        self._linelist_cache = linelist




    def make_linelist(self, starname: str):
        """Run ARES, merge, build EW tables, compute EW errors for <4-line elements."""
        # Prepare ARES inputs in CWD
        self._write_cdo_from_linelist()
        self._write_mine_opt(starname)

        # Run ARES (must be installed and available in PATH)
        os.system('ARES')

        # Load inputs/outputs — reuse cached linelist to guarantee identical row count
        linelist = getattr(self, "_linelist_cache", None)
        if linelist is None:
            linelist = pd.read_csv(
                self.linelist_path,
                skiprows=2,
                names=['WL', 'EP', 'loggf', 'element', 'num', 'damp', 'EWsun'],
                sep=r'\s+',
                engine='python',
                dtype={'WL': float},
            ).sort_values(['WL'])

        aresout = pd.read_csv(
            "aresout.dat",
            usecols=(0,2,3,4,5),
            names=['wave', 'depth', 'FWHM', 'EW', 'EWerr'],
            sep=r'\s+',
            engine='python',
            dtype={'wave': float},
        ).sort_values(['wave'])

        # ----- robust merge: match on wavelength rounded to 2 decimals on BOTH sides -----
        linelist = linelist.copy()
        aresout = aresout.copy()
        linelist['WL_2d'] = linelist['WL'].round(2)
        aresout['wave_2d'] = aresout['wave'].round(2)

        data_output = pd.merge(
            linelist,
            aresout,
            left_on='WL_2d',
            right_on='wave_2d',
            how='inner',
        )

        # (optional) keep original columns clean
        data_output = data_output.drop(columns=['WL_2d', 'wave_2d'])

        # Filter bad depths (as before); if you want to keep everything, comment this out
        data_output = data_output[(data_output.depth >= 0) & (data_output.depth < 1.0)]

        data_output_path = self.dir_out / f"{starname}_aresout.dat"
        data_output.to_csv(data_output_path, index=False, sep='\t')

        # EW MOOG input
        log_path = Path("logARES.txt")
        snr = self._read_snr_from_log(log_path)
        hdr = f'{starname} - S/N: {snr}'
        ew_moog_input = data_output.loc[:, ['WL', 'num', 'EP', 'loggf', "damp", 'EW']].sort_values(['num'])
        fileout = Path(f"aresout_{starname}.dat")
        np.savetxt(fileout.as_posix(), ew_moog_input.values,
                   fmt=('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%9.2f', '%19.1f'),
                   header=f' {hdr}')

        # Select elements with <4 lines
        element_list_tmp = data_output.groupby('num')['num'].count()
        element_list = pd.DataFrame([element_list_tmp.index, element_list_tmp.values],
                                    index=['num', 'nlines']).transpose()
        elements_with_lt_4_lines = element_list[element_list.nlines < 4]
        lines_lt4 = pd.merge(left=data_output, right=elements_with_lt_4_lines, on='num')

        # Cayrel error
        lines_lt4['ew_plus_error'] = np.round(
            lines_lt4['EW'] + 1000*2.45*((lines_lt4['FWHM']*0.035/2.3548)**0.5)/snr + 6*1000*lines_lt4['FWHM']/snr/2.3548, 1
        )
        lines_lt4['ew_minus_error'] = np.round(
            lines_lt4['EW'] - (1000*2.45*((lines_lt4['FWHM']*0.035/2.3548)**0.5)/snr + 6*1000*lines_lt4['FWHM']/snr/2.3548), 1
        )

        # EW ± error tables
        ew_plus = lines_lt4.loc[:, ['WL', 'num', 'EP', 'loggf', "damp", 'ew_plus_error']].rename(columns={"ew_plus_error": "ew_error"})
        ew_minus = lines_lt4.loc[:, ['WL', 'num', 'EP', 'loggf', "damp", 'ew_minus_error']].rename(columns={"ew_minus_error": "ew_error"})
        ew_pm = pd.concat([ew_plus, ew_minus], ignore_index=True)
        ew_pm = ew_pm[ew_pm['ew_error'] > 0].sort_values(['num'])

        fileout_err = Path(f"aresout_{starname}_ew_err.dat")
        np.savetxt(fileout_err.as_posix(), ew_pm.values,
                   fmt=('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%9.2f', '%19.1f'),
                   header=f' {hdr}')

        # Move logs and EWs to outputs
        (self.dir_log / f"{starname}.log").write_text(log_path.read_text())
        data_output_path = self.dir_out / f"{starname}_aresout.dat"  # already written
        (self.dir_ew / f"{starname}_aresout.dat").write_bytes(fileout.read_bytes())
        (self.dir_ew_err / f"{starname}_aresout_ew_err.dat").write_bytes(fileout_err.read_bytes())

        # cleanup CWD temp files
        for tmp in ["cdo.dat", "mine.opt", "aresout.dat", fileout.name, fileout_err.name, "logARES.txt"]:
            try:
                Path(tmp).unlink()
            except FileNotFoundError:
                pass

    def ew_measurements(self):
        """Run ARES for all stars in the input DataFrame and materialize outputs/* directories."""
        for i in range(len(self.input_df)):
            star = str(self.input_df.star[i])
            print(star)
            self.make_linelist(star)

        # Create/refresh outputs/ares_ew_corr by copying the *contents* of ares_ew
        import shutil
        self.dir_corr.mkdir(parents=True, exist_ok=True)
        for item in self.dir_ew.iterdir():
            dst = self.dir_corr / item.name
            if item.is_dir():
                shutil.copytree(item, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(item, dst)
        print("Synced outputs/ares_ew -> outputs/ares_ew_corr (contents only).")
        print("Done")

