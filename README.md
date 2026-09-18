# MultiFischer

[![CI](https://github.com/HamishGTrowell/MultiFischer/actions/workflows/ci.yml/badge.svg)](https://github.com/HamishGTrowell/MultiFischer/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

MultiFischer is a command-line program for estimating photoswitch photostationary state (PSS) compositions from UV-Vis spectra. It extends the two-irradiation [Fischer method](https://doi.org/10.1021/j100870a063) by applying every unique irradiation wavelength pair, aggregating the resulting PSS estimates, and filtering nonphysical and anomalous results.

> [!NOTE]
> The multipair sensitivity analysis has a known issue. An update is planned for a future release.

## What MultiFischer does

- calculates standard Fischer PSS estimates using the quantum yield ratio invariance assumption (`X = 1`);
- uses all unique irradiation pairs and imputes each pair's result at every observed irradiation wavelength;
- applies pair-level PSS range and metastable spectrum filters followed by an irradiation-specific empirical HDI filter;
- reports means, sample standard deviations, accepted-value counts, and full intermediate results;
- produces CSV results tables, diagnostic logs, and plots; and
- supports a single-pair sensitivity mode for inspecting two PSS curves across the selected quantum-yield-ratio range.

## Installation with Conda

Clone the repository, create the environment, and activate it:

```bash
git clone https://github.com/HamishGTrowell/MultiFischer.git
cd MultiFischer
conda env create -f environment.yml
conda activate multifischer
```

The environment installs Python 3.12, the scientific Python dependencies, and MultiFischer itself in editable mode. For a non-editable installation, replace `-e .` in `environment.yml` with `.` before creating the environment.

## Quick start

Run the included example with the default multipair analysis:

```bash
multifischer --input example.csv --out example-output
```

Run one Fischer pair across the default quantum-yield-ratio range:

```bash
multifischer --input example.csv --out example-singlepair --singlepair 365 385
```

All options are listed by:

```bash
multifischer --help
```

## Input format

Input is a comma-separated file with:

1. a wavelength column in nm, headed `Wavelength` (common variants such as `Wavelength (nm)` and `WL` are accepted);
2. the spectrum of the dark or stable state, headed `Dark`; and
3. at least two columns of PSS UV-Vis absorbance data whose headings are their irradiation wavelengths, for example, `365` or `385 nm`.

```csv
Wavelength,Dark,365,385,400
650.0,0.0031,0.0030,0.0032,0.0030
649.0,0.0031,0.0030,0.0032,0.0030
...
```

The supplied [`example.csv`](example.csv) is a complete working example. See the [user guide](docs/user-guide.md) for validation rules, additional options, and output-file descriptions.

## Scientific basis

The original Fischer method estimates conversion in a reversible two-state photochemical system from a known dark-state spectrum and PSS spectra obtained using two irradiation wavelengths. Its central assumption is that the forward-to-reverse quantum-yield ratio is invariant between the irradiation wavelengths.

MultiFischer treats `X = 1` as the standard invariant-ratio result and applies the calculation over every unique wavelength pair. Single-pair mode evaluates a selected irradiation pair over a configurable range of `X` values.

> [!NOTE]
> The multipair sensitivity analysis has a known issue. An update is planned for a future release.

The exact equations, pair aggregation, filters, summary statistics, limitations, and interpretation are documented in the [scientific guide](docs/scientific-method.md).

## Reproducibility

Every run writes `run_settings.json`, recording the command, normalized settings, timestamp, Python version, operating system, executable, working directory, and package versions. The accompanying `multifischer.log` records the calculation and filter stages. Retain both files with exported CSV results.

## Citation

Use the repository's [`CITATION.cff`](CITATION.cff) through GitHub's “Cite this repository” panel. Scientific work should also cite the original method:

> E. Fischer, “Calculation of photostationary states in systems A ⇄ B when only A is known,” *J. Phys. Chem.* **71** (1967), 3704–3706. https://doi.org/10.1021/j100870a063

## Contributing

Bug reports, proposed tests, and documentation improvements are welcome. See [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License

MultiFischer is released under the permissive [MIT License](LICENSE). This permits reuse, modification, and redistribution with preservation of the copyright and license notice.
