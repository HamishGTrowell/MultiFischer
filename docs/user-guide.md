# xFischer user guide

## 1. Installation

### Recommended Conda environment

From a clone of the repository:

```bash
conda env create -f environment.yml
conda activate xfischer
xfischer --help
```

To construct the same type of environment manually:

```bash
conda create -n xfischer -c conda-forge python=3.12 numpy pandas matplotlib scipy seaborn pip
conda activate xfischer
conda env config vars set MPLBACKEND=Agg
conda deactivate
conda activate xfischer
python -m pip install -e .
```

The editable installation is convenient when working from a clone. Use `python -m pip install .` for a regular installation. `MPLBACKEND=Agg` selects a non-interactive backend suitable for xFischer's file-based plots and for headless systems. The module form `python -m xfischer` is equivalent to the installed `xfischer` command.

## 2. Preparing input data

xFischer accepts one CSV file per analysis. Rows are spectral sampling points and columns are spectra.

Required columns:

- **Wavelength:** wavelength in nm. Values must be numeric, finite, unique, and contain at least two rows. Row order does not matter.
- **Dark:** absorbance of the dark-adapted or stable state. Values must be numeric and finite.
- **Irradiated PSS spectra:** at least two absorbance columns. Each heading must be a numeric irradiation wavelength, optionally followed by `nm`, such as `365`, `365.0`, or `365 nm`. Values must be numeric and finite.

Headings are case-insensitive and surrounding whitespace is ignored. Accepted wavelength aliases include `wavelength`, `wavelength (nm)`, `wavelength/nm`, `wl`, and `nm`. Accepted dark-state aliases include `dark`, `ambient`, `stable`, `pure E`, `pure Z`, `non-irradiated`, and `initial`.

Irradiation wavelengths must fall inside the measured spectral range. The requested dark-maximum, metastate-filter, and plotting ranges must each overlap the measured wavelengths. xFischer uses the measured row nearest to each irradiation wavelength; interpolation is not applied at this extraction step.

Example header:

```csv
Wavelength,Dark,365,385,400,420,433,450,505,530,590
```

## 3. Running an analysis

### Multipair analysis

```bash
xfischer --input example.csv --out example-output
```

All unique irradiation pairs are evaluated at `X = 1`, filtered, summarized, and then reevaluated at 200 evenly spaced `X` values over the sensitivity range.

If `--out` is omitted, the output directory is named after the input file stem. If the directory already exists, xFischer asks whether to reuse it or select another path. Existing files with standard output names may be replaced when an existing directory is reused.

### Single-pair analysis

```bash
xfischer \
  --input example.csv \
  --out example-singlepair \
  --singlepair 365 385 \
  --qyratio-range 0.5 2.0
```

Single-pair mode calculates the PSS at both selected irradiation wavelengths over the ratio range. Multipair filters and summary plots are not produced.

## 4. Command-line options

| Option | Default | Meaning |
|---|---:|---|
| `--input INPUT.CSV` | required | Input UV-Vis CSV file. |
| `--out OUT/` | input stem | Output directory. |
| `--singlepair IRR1 IRR2` | off | Analyze only one irradiation pair. |
| `--qyratio-range MIN_X MAX_X` | `0.5 2.0` | Positive range of quantum-yield-ratio values for sensitivity analysis. |
| `--darkmax-range MIN_NM MAX_NM` | `275 800` | Band in which the maximum of the dark spectrum is located. |
| `--meta-range MIN_NM MAX_NM` | `275 600` | Spectral band used by the metastable-spectrum physicality filter. |
| `--meta-min MIN_ABS` | `-0.05` | Minimum permitted extrapolated metastable absorbance in the filter band. |
| `--pss-range MIN_PSS MAX_PSS` | `-0.05 1.05` | Inclusive accepted range for fractional PSS values. |
| `--hdi HDI` | `0.9` | Fraction of candidate values in the narrowest empirical HDI core. |
| `--hdi-sigma HDI_SIGMA` | `5` | Sample-standard-deviation multiple used to extend the HDI core bounds. |
| `--uvvis-range MIN_NM MAX_NM` | `200 600` | Horizontal range of the UV-Vis plot. |
| `--violin-irrs IRR [...]` | automatic | Irradiations to show in the violin plot; use `all` for every irradiation. |
| `--sensitivity` | off | Add sensitivity ranges to the final UV-Vis and PSS plots. |
| `--metastate LABEL` | `Z` | Label used for the metastable isomer in figures. |

PSS values are represented internally and in CSV files as fractions. Plot axes generally display percentages.

## 5. Multipair outputs

| File | Contents |
|---|---|
| `run_settings.json` | Command, normalized settings, timestamp, platform, Python executable, and dependency versions. |
| `xfischer.log` | Human-readable run log, diagnostics, filter counts, and summary values. |
| `pss_matrix.csv` | Every pair/observed-irradiation PSS at `X = 1`, including filter decisions. |
| `metastate_spectra.csv` | Pair-specific extrapolated metastable spectra in wide form. |
| `pss_accepted_vals.csv` | Rows from the PSS matrix passing every filter. |
| `hdi_results.csv` | Candidate counts, HDI core statistics, extended bounds, and pass/fail counts by irradiation. |
| `pss_final_summary.csv` | Final mean PSS, sample standard deviation, and accepted count by irradiation at `X = 1`. |
| `pss_sensitivity_summary.csv` | Independently filtered summary at every sampled `X` and irradiation wavelength. |
| `pair_sensitivity_matrix.csv` | Directed pairwise PSS range across `X`, in percentage points. |
| `pair_error_matrix.csv` | Directed direct-pair deviation from the final mean at `X = 1`, in percentage points. |
| `uv-vis.png` | Dark, irradiated, and representative metastable spectra. |
| `pss_vs_irradiation.png` | Final PSS summary against irradiation wavelength. |
| `pss_distribution_*.png` | PSS-distribution violin plot for the selected irradiation(s). |
| `multipair_qy_sensitivity.png` | Filtered mean PSS against `X`. |
| `pair_sensitivity_heatmap.png` | Directed pair sensitivity heatmap. |
| `pair_error_heatmap.png` | Directed pair error heatmap. |
| `error_sensitivity_scatter.png` | Pair error versus ratio sensitivity. |

With `--sensitivity`, the UV-Vis and PSS summary figures include envelopes derived from the sampled quantum-yield-ratio results.

## 6. Single-pair outputs

For `--singlepair 365 385`, the run produces:

- `365_385nm_qy_ratio_sensitivity.csv`;
- `365_385nm_qy_ratio_sensitivity.png`;
- `run_settings.json`; and
- `xfischer.log`.

## 7. Interpreting diagnostics

Review `pss_matrix.csv` and `hdi_results.csv` rather than relying only on final plots. A small accepted count, rejection of most pairs, strong variation across `X`, or a direct pair with high error and high sensitivity indicates that the reported mean is fragile for that irradiation wavelength.

If no finite values pass every filter at `X = 1`, the run stops with an explanatory error. Inspect input spectra, confirm that each spectrum represents a true PSS, and review the filter settings before relaxing thresholds.

## 8. Reproducible reporting checklist

Archive the input CSV, xFischer version, complete output directory, and `run_settings.json`. Report the irradiation wavelengths, all non-default ranges and filters, accepted counts, final mean and sample standard deviation, and the tested `X` range. State explicitly that the sensitivity range is an assumption analysis rather than an experimentally determined quantum-yield interval.
