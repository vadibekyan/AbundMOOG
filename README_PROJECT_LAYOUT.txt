
Project layout (proposed & applied in this session)
==================================================

inputs/
  - spectra/                             # (optional) raw spectra
  - input_param_error.rdb                # star list + stellar params (move here)
  - elements_formatting.rdb              # element list + solar values
  - elements_to_plot.rdb                 # which elements to include in plots
  - abundances_harps1111_o_c.rdb         # HARPS comparison sample

src/
  - myspec.py                            # pipeline entry-point (run with PYTHONPATH=src)
  - ares_ew.py
  - atmospheres_kurucz.py
  - abundances_3sigma_outlier.py
  - abund_results_formatting.py
  - elfe_plot.py

outputs/
  - ares_ew_corr/                        # ARES-corrected EWs (and outputs/ares_ew_corr/err_ew/)
  - atmospheres/                         # atmosphere .atm files
  - moog_abundances/
      moog_outab/
      moog_std_out/
      moog_outab_formatted/
      moog_abstar/
  - analysis/
      all_abundances_ordered_element.rdb
      all_abundances_ordered_star.rdb
  - plots/
      elfe_feh.pdf, elfe_teff.pdf, elfe_logg.pdf, HR.pdf

How to run
----------
# from project root (where inputs/ and src/ sit)
export PYTHONPATH=src
python src/myspec.py

Notes
-----
- Paths were updated in the modules to point to inputs/* and outputs/* by default.
- You can still override directories when instantiating the classes if needed.


Run without PYTHONPATH (recommended)
-----------------------------------
# From project root (where src/ and inputs/ live)
python -m src.myspec --config config.json --target-name "TOI-2128_FIES2020"

Notes:
- `config.json` lets you change input/output files and plotting target without editing code.
