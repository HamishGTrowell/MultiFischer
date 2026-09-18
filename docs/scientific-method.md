# Scientific method and interpretation

## Scope

MultiFischer analyzes a reversible two-state photoswitch system using a known dark/stable-state UV-Vis spectrum and spectra measured at photostationary state (PSS) under at least two irradiation wavelengths. It implements the Fischer calculation, extends it over all unique irradiation pairs, applies physical and empirical filters, and provides single-pair sensitivity analysis.

> [!NOTE]
> The multipair sensitivity analysis has a known issue. An update is planned for a future release.

The original method is described by E. Fischer, [“Calculation of photostationary states in systems A ⇄ B when only A is known,” *J. Phys. Chem.* **71** (1967), 3704–3706](https://doi.org/10.1021/j100870a063).

## Notation and pair calculation

For an ordered irradiation pair `(λ₁, λ₂)`, let:

- `A_D(λ)` be the dark-state absorbance at wavelength `λ`;
- `A_i(λ)` be the spectrum measured at PSS under irradiation `λᵢ`;
- `λ_max` be the measured wavelength with maximum `A_D` inside `--darkmax-range`;
- `d₁ = A₁(λ₁) − A_D(λ₁)` and `d₂ = A₂(λ₂) − A_D(λ₂)`; and
- `n = [A₂(λ_max) − A_D(λ_max)] / [A₁(λ_max) − A_D(λ_max)]`.

The experimental row nearest each irradiation wavelength supplies the absorbance at `λ₁` or `λ₂`.

At the standard Fischer condition `X = 1`, the fractional PSS conversion for the first irradiation is implemented as

```text
PSS₁ = [d₂/A_D(λ₂) − d₁/A_D(λ₁)]
       -------------------------------------------------------
       [1 + d₂/A_D(λ₂) − n(1 + d₁/A_D(λ₁))].
```

For the second irradiation, pair order is reversed.

## Single-pair sensitivity calculation

`X` is a dimensionless perturbation of the effective forward-to-reverse quantum-yield-ratio relationship between the two irradiation wavelengths. The standard invariance result is `X = 1`. In single-pair sensitivity analysis, MultiFischer solves the following equation for the first irradiation at general `X`:

```text
a PSS₁² + b PSS₁ + c = 0
```

with

```text
a = nX − n
b = 1 − n d₁/A_D(λ₁) − nX + X d₂/A_D(λ₂)
c = d₁/A_D(λ₁) − X d₂/A_D(λ₂).
```

Of the two quadratic roots, the implementation selects the branch continuous at `X = 1`, using the sign of the limiting linear coefficient. Degenerate linear cases are solved directly; undefined divisions and negative discriminants produce missing values. For the second irradiation, pair order is reversed and `1/X` is used.

The default `--qyratio-range 0.5 2.0` is sampled at 200 evenly spaced points including both endpoints. This range is a sensitivity scenario, not a confidence interval and not an estimate of quantum yields.

## Multipair extension

For `m` irradiation spectra, MultiFischer constructs all `m(m−1)/2` unique pairs. Each pair provides direct PSS estimates at its two irradiation wavelengths. At every remaining observed irradiation wavelength `λ₃`, the pair estimate is linearly interpolated or extrapolated along absorbance at `λ_max`:

```text
PSS₃ = PSS₁ + [A₃(λ_max) − A₁(λ_max)]
                  × (PSS₂ − PSS₁)
                  / [A₂(λ_max) − A₁(λ_max)].
```

Thus the unfiltered PSS matrix contains one row for every pair and observed irradiation combination. Equal pair absorbances at `λ_max` make third-wavelength imputation undefined.

## Extrapolated metastable-state spectra

For a pair result `PSSᵢ`, each observed spectrum independently implies a pure metastable-state spectrum:

```text
A_meta,i(λ) = A_D(λ) + [A_i(λ) − A_D(λ)] / PSSᵢ.
```

The pair-level metastable spectrum used for filtering is the pointwise arithmetic mean of the two independently extrapolated spectra. Zero or undefined PSS denominators make that pair's metastable spectrum unavailable.

For the final UV-Vis figure, MultiFischer selects the irradiation wavelength with the highest accepted mean PSS. It extrapolates a representative metastable spectrum using the mean PSS and an absorbance envelope using mean PSS ± one sample standard deviation.

## Filters

For the standard multipair result at `X = 1`, filters are applied sequentially.

### 1. Pair-level PSS range

Every PSS row belonging to a pair must be finite and within the inclusive `--pss-range` (default `−0.05` to `1.05`). If any observed irradiation estimate fails, the entire pair fails.

### 2. Pair-level metastable-spectrum filter

For each pair, MultiFischer finds the minimum extrapolated metastable absorbance within `--meta-range` (default 275–600 nm). The pair passes when that minimum is at least `--meta-min` (default −0.05 AU). This rejects strongly negative, physically implausible extrapolations while allowing modest baseline noise.

### 3. Irradiation-specific empirical HDI filter

For each observed irradiation wavelength separately, candidate PSS values are those that passed both earlier filters. Values are sorted and the narrowest contiguous window containing at least `ceil(hdi × N)` values is selected as the empirical HDI core. With core mean `μ_core` and sample standard deviation `s_core`, the acceptance interval is

```text
[μ_core − hdi_sigma × s_core, μ_core + hdi_sigma × s_core].
```

Defaults are `hdi = 0.9` and `hdi_sigma = 5`. The HDI here is an empirical shortest interval used for outlier filtering; it is not a Bayesian credible interval. A one-value core has undefined sample standard deviation, so no value passes its HDI bounds.

## Final statistics

At each irradiation wavelength, the reported mean, sample standard deviation (`ddof = 1`), and count are calculated from all rows passing all three filters. The HDI core defines filter bounds only and is not itself the reported sample.

At `X = 1`, these are the standard final PSS results.

## Multipair sensitivity analysis

> [!NOTE]
> The multipair sensitivity analysis has a known issue. An update is planned for a future release.

## Pair diagnostics

The directed pair error matrix reports, in percentage points, the absolute difference between a direct pair's `X = 1` PSS and the final accepted multipair mean at the same target irradiation wavelength.

## Assumptions and limitations

Interpretation requires all of the following to be reasonable for the experiment:

- the system is adequately described by two interconverting absorbing states;
- the dark spectrum represents the pure or known stable state used by the calculation;
- every irradiated spectrum was measured after reaching PSS;
- spectra share concentration, optical path length, solvent, temperature, baseline treatment, and instrumental response;
- Beer-Lambert additivity is applicable over the measured range;
- thermal interconversion, photodegradation, side reactions, aggregation, and concentration drift do not materially distort the PSS spectra; and
- the effective forward/reverse quantum-yield-ratio relationship is sufficiently invariant for the standard Fischer result.

Multiple pair estimates are not independent experimental replicates because they reuse the same spectra. Their standard deviation characterizes dispersion among pair-derived estimates after filtering; it is not automatically an experimental standard error or a complete uncertainty budget.

## Recommended reporting

Report the MultiFischer version, input irradiation wavelengths, `λ_max` search band and selected wavelength, all filter settings, number of accepted estimates per irradiation, and mean and sample standard deviation at `X = 1`. Retain the input spectra, `run_settings.json`, `multifischer.log`, full PSS matrix, accepted rows, and HDI diagnostics with the publication archive.
