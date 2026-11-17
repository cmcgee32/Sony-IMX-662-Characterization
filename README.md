# Sony STARVIS CMOS X-ray characterization (IMX662, IMX226, IMX290)

MATLAB (and Python table) code used by Carissma McGee (MIT Kavli Institute, 2026) to characterize
Sony STARVIS CMOS sensors for soft X-ray astronomy, following the lab method of Lupo (2025, MIT S.M. thesis)
and compared with Schneider et al. (2024, JATIS 10(3) 038001).

## What is measured
| Test | Scripts |
|---|---|
| Hot pixels, X-ray spectra, energy resolution (Fe-55, Po-210) | `imx662_xray_events.m`, `imx662_plot_xray.m`, `imx226_xray_from_reduced.m` |
| Light curves and count-rate heat maps | `imx662_plot_xray.m` |
| Readout noise | `imx662_readnoise.m`, `imx226_readnoise.m` |
| Gain and pixel-to-pixel gain variation | `imx662_xray_events.m`, `imx290_ecc_analysis.m` |
| X-ray linearity and energy-resolution budget | `imx662_xray_linearity.m`, `imx226_tables.m` |
| Optical linearity and photon transfer | `imx662_optical_linearity.m`, `imx226_optical_linearity.m` |
| Dark current | `imx662_dark_current.m`, `imx226_dark_current.m` |
| IMX226 polarimetry-beamline event rates | `imx226_beamline.m` |
| IMX290 per-pixel gain calibration (ECC) | `imx290_ecc_analysis.m`, `imx290_plot.m` |
| Summary tables (booktabs LaTeX/PDF/PNG) | `imx662_tables.m`, `make_tables.py`, `make_tables_new.py` |

## How to run
1. Put the lab data zips (e.g. `IMX662_Fe55_SpectrumData_20250625.zip`) in the folder set in `imx662_cfg.m` (`c.D`).
2. In MATLAB: `stages = [3 4 1 2 5 6]; imx662_run_all` for the IMX662, and `imx226_run` for the IMX226.
3. Figures (PNG + FIG), CSV tables and `.mat` results are written to `IMX662_results/` and `IMX226_results/`.
4. `python3 make_tables.py <results> <out>` builds the academic (booktabs) tables.

## Key results
| | IMX662 | IMX226 | IMX290 |
|---|---|---|---|
| Conversion gain [ADU/eV] | 0.2652 | 0.3646 | 0.2941 |
| Readout noise, mode [e- rms] | 3.53 | 2.30 | 1.9 (Schneider) |
| Mn K-alpha FWHM, single-pixel [eV] | 133.3 | 140.9 | 134.2 (122.7 with per-pixel gain) |
| Pixel-to-pixel gain variation [%] | 0.24 | 0.38 | 0.40 |
| Dark current, lab method [e-/pixel/sec] | 0.027 | 0.438 | - |
