# xFischer

[![CI](https://github.com/HamishGTrowell/xFischer/actions/workflows/ci.yml/badge.svg)](https://github.com/HamishGTrowell/xFischer/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

xFischer is a command-line program for estimating photoswitch photostationary-state (PSS) compositions from UV-Vis spectra. It extends the two-irradiation [Fischer method](https://doi.org/10.1021/j100870a063) by applying every unique irradiation-wavelength pair, aggregating the resulting PSS estimates, filtering nonphysical and anomalous results, and testing sensitivity to the quantum-yield-ratio invariance assumption.

## What xFischer does

- calculates standard Fischer PSS estimates at a quantum-yield ratio of `X = 1`;
- uses all unique irradiation pairs and imputes each pair's result at every observed irradiation wavelength;
- applies pair-level PSS-range and metastable-spectrum filters followed by an irradiation-specific empirical HDI filter;
- reports means, sample standard deviations, accepted-value counts, and full intermediate results;
- repeats the filtered multipair calculation over 200 quantum-yield-ratio values to quantify assumption sensitivity;
- produces publication-ready CSV tables, diagnostic logs, and plots; and
- supports a single-pair mode for inspecting two PSS curves across the selected ratio range.

## Installation with Conda

Clone the repository, create the environment, and activate it:

```bash
git clone https://github.com/HamishGTrowell/xFischer.git
cd xFischer
conda env create -f environment.yml
conda activate xfischer
```

The environment installs Python 3.12, the scientific Python dependencies, and xFischer itself in editable mode. It also selects Matplotlib's non-interactive `Agg` backend because xFischer saves figures to files rather than opening GUI windows. For a non-editable installation, replace `-e .` in `environment.yml` with `.` before creating the environment.

## Quick start

Run the included example with the default multipair analysis:

```bash
xfischer --input example.csv --out example-output
```

Include sensitivity envelopes on the final UV-Vis and PSS plots:

```bash
xfischer --input example.csv --out example-output --sensitivity
```

Run one Fischer pair across the default quantum-yield-ratio range:

```bash
xfischer --input example.csv --out example-singlepair --singlepair 365 385
```

All options are listed by:

```bash
xfischer --help
```

## Input format

Input is a comma-separated file with:

1. a wavelength column in nm, headed `Wavelength` (common variants such as `Wavelength (nm)` and `WL` are accepted);
2. the spectrum of the dark or stable state, headed `Dark` (several common aliases are accepted); and
3. at least two PSS spectra whose headings are their irradiation wavelengths, for example `365`, `385`, or `400 nm`.

```csv
Wavelength,Dark,365,385,400
650.0,0.0031,0.0030,0.0032,0.0030
649.0,0.0031,0.0030,0.0032,0.0030
...
```

The supplied [`example.csv`](example.csv) is a complete working example. See the [user guide](docs/user-guide.md) for validation rules, every option, and output-file descriptions.

## Scientific basis

The original Fischer method estimates conversion in a reversible two-state photochemical system from a known dark-state spectrum and PSS spectra obtained using two irradiation wavelengths. Its central assumption is that the forward-to-reverse quantum-yield ratio is invariant between the irradiation wavelengths.

xFischer treats `X = 1` as the standard invariant-ratio result, applies the calculation over every unique wavelength pair, and repeats the full filtered analysis over a configurable range of `X` values. This sensitivity analysis shows how strongly the reported PSS depends on deviations from the invariance assumption; it does not measure the quantum yields or assign a probability distribution to `X`.

The exact equations, pair aggregation, filters, summary statistics, limitations, and interpretation are documented in the [scientific guide](docs/scientific-method.md).

## Reproducibility

Every run writes `run_settings.json`, recording the command, normalized settings, timestamp, Python version, operating system, executable, working directory, and package versions. The accompanying `xfischer.log` records the calculation and filter stages. Retain both files with exported CSV results.

## Citation

Use the repository's [`CITATION.cff`](CITATION.cff) through GitHub's “Cite this repository” panel. Scientific work should also cite the original method:

> E. Fischer, “Calculation of photostationary states in systems A ⇄ B when only A is known,” *J. Phys. Chem.* **71** (1967), 3704–3706. https://doi.org/10.1021/j100870a063

## Contributing

Bug reports, proposed tests, and documentation improvements are welcome. See [`CONTRIBUTING.md`](CONTRIBUTING.md). Changes that affect scientific calculations should include a mathematical rationale and regression tests.

## License

xFischer is released under the permissive [MIT License](LICENSE). This permits reuse, modification, and redistribution with preservation of the copyright and license notice.
