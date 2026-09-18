"""Core PSS calculations for MultiFischer."""

from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from itertools import combinations

from . import utils
from . import outlog
from . import columns as cols

N_QY_POINTS = 200


@dataclass(frozen=True)
class AbsorbanceData:
    """Store absorbance values required for one Fischer-pair calculation.

    Parameters
    ----------
    dark_irr1, dark_irr2 : float
        Dark-spectrum absorbance at irradiation wavelengths ``irr1`` and
        ``irr2``.
    dark_max : float
        Dark-spectrum absorbance at its selected maximum.
    irr1_max, irr2_max : float
        Irradiated-spectrum absorbance at the selected dark-spectrum maximum.
    irr1_at_irr1, irr2_at_irr2 : float
        Absorbance of each irradiated spectrum at its own irradiation
        wavelength.
    """
    dark_irr1: float
    dark_irr2: float
    dark_max: float
    irr1_max: float
    irr2_max: float
    irr1_at_irr1: float
    irr2_at_irr2: float

    @property
    def d1(self) -> float:
        """Return the irradiation-induced absorbance change at ``irr1``."""
        return self.irr1_at_irr1 - self.dark_irr1

    @property
    def d2(self) -> float:
        """Return the irradiation-induced absorbance change at ``irr2``."""
        return self.irr2_at_irr2 - self.dark_irr2

    def reversed(self) -> 'AbsorbanceData':
        """Return a copy with the two irradiation spectra exchanged."""

        return AbsorbanceData(
            dark_irr1=self.dark_irr2,
            dark_irr2=self.dark_irr1,
            dark_max=self.dark_max,
            irr1_max=self.irr2_max,
            irr2_max=self.irr1_max,
            irr1_at_irr1=self.irr2_at_irr2,
            irr2_at_irr2=self.irr1_at_irr1,
        )


@dataclass
class FischerPairAnalysis:
    """Perform Fischer and MultiFischer calculations for one irradiation pair.

    Absorbance values needed by the Fischer equations are extracted during
    initialisation. Calculated PSS values are cached by quantum-yield ratio so
    that repeated imputation, metastate extrapolation, and sensitivity
    calculations do not repeat the same pair calculation.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Validated and normalised UV–Vis data.
    irr1, irr2 : str
        Irradiation wavelength labels defining the Fischer pair.
    darkmax_range : tuple[float, float]
        Wavelength range, in nm, used to locate the dark-spectrum absorbance
        maximum.
    x_array : numpy.ndarray
        Quantum-yield-ratio values used for sensitivity analysis.

    Attributes
    ----------
    abs_data : AbsorbanceData
        Absorbance values extracted for the Fischer calculation.
    idx_darkmax : int
        DataFrame index of the selected dark-spectrum absorbance maximum.
    """
    data_df: pd.DataFrame
    irr1: str
    irr2: str
    darkmax_range: tuple[float, float]
    x_array: np.ndarray

    abs_data: AbsorbanceData = field(init=False)
    idx_darkmax: int = field(init=False)

    _pss_cache: dict[float, tuple[float, float]] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )

    def __post_init__(self):
        """Locate the dark-state maximum and extract pair absorbance values."""
        self.idx_darkmax, _ = utils.dark_lambda_max(
            self.data_df,
            self.darkmax_range,
        )
        
        self.abs_data = extract_abs_data(
            self.data_df,
            self.irr1,
            self.irr2,
            self.idx_darkmax,
        )
    
    def pss_values(self, X=1.0):
        """Calculate and cache PSS values for both pair wavelengths.

        Parameters
        ----------
        X : float, optional
            Quantum-yield ratio. The default of 1 gives the standard Fischer
            result.

        Returns
        -------
        tuple[float, float]
            PSS values at irradiation wavelengths ``irr1`` and ``irr2``.
        """

        cache_key = float(X)

        if cache_key not in self._pss_cache:
            self._pss_cache[cache_key] = (
                multifischer(self.abs_data, cache_key, target='irr1'),
                multifischer(self.abs_data, cache_key, target='irr2'),
            )

        return self._pss_cache[cache_key]

    def irr1_pss(self, X=1.0):
        """Return the calculated PSS at irradiation wavelength ``irr1``."""
        return self.pss_values(X)[0]

    def irr2_pss(self, X=1.0):
        """Return the calculated PSS at irradiation wavelength ``irr2``."""
        return self.pss_values(X)[1]

    def impute_irr3_pss(self, irr3, X=1.0):
        """Impute the PSS at a third observed irradiation wavelength.

        Linear interpolation or extrapolation is performed from the pair PSS values
        and the irradiated-spectrum absorbances at the selected dark-spectrum
        maximum.

        Parameters
        ----------
        irr3 : str
            Observed irradiation wavelength label to impute.
        X : float, optional
            Quantum-yield ratio.

        Returns
        -------
        float
            Imputed PSS at ``irr3``, or ``numpy.nan`` when both pair spectra
            have the same absorbance at the selected dark-state maximum.
        """
        a_irr3_max = self.data_df.at[self.idx_darkmax, irr3]

        pss1, pss2 = self.pss_values(X)

        absorbance_span = (
            self.abs_data.irr2_max
            - self.abs_data.irr1_max
        )

        # Equal pair absorbances at the dark-state maximum provide no axis
        # along which a third-wavelength PSS can be interpolated
        if absorbance_span == 0.0:
            return np.nan

        pss3 = pss1 + (
            (a_irr3_max - self.abs_data.irr1_max)
            * (pss2 - pss1)
            / absorbance_span
        )

        return pss3

    def extrapolate_meta_spectrum(self, X=1.0):
        """Extrapolate the metastable-state absorbance spectrum for the pair.

        A metastate spectrum is independently extrapolated from each irradiated
        spectrum. At every spectral wavelength, the arithmetic mean of the two
        extrapolated absorbances is retained.

        Parameters
        ----------
        X : float, optional
            Quantum-yield ratio.

        Returns
        -------
        pandas.DataFrame
            Long-form metastate spectrum containing the Fischer-pair identifiers,
            spectral wavelength, and extrapolated absorbance. The extrapolated
            absorbance is missing when either pair PSS is undefined or zero.
        """
        pss1, pss2 = self.pss_values(X)

        # Undefined or zero PSS values cannot be used as extrapolation
        # denominators. Mark this pair unavailable so the pair-level filters
        # reject it without producing divide-by-zero warnings or infinities
        if (
            not np.isfinite(pss1)
            or not np.isfinite(pss2)
            or pss1 == 0.0
            or pss2 == 0.0
        ):
            meta_abs_mean = np.full(
                len(self.data_df),
                np.nan,
                dtype=float,
            )
        else:
            meta_abs_irr1 = (
                self.data_df[cols.DARK]
                + (self.data_df[self.irr1] - self.data_df[cols.DARK]) / pss1
            )
            meta_abs_irr2 = (
                self.data_df[cols.DARK]
                + (self.data_df[self.irr2] - self.data_df[cols.DARK]) / pss2
            )

            meta_abs_mean = np.mean(
                (
                    meta_abs_irr1.to_numpy(dtype=float),
                    meta_abs_irr2.to_numpy(dtype=float),
                ),
                axis=0,
            )

        return pd.DataFrame(
            {
                cols.PAIR_IRR1: self.irr1,
                cols.PAIR_IRR2: self.irr2,
                cols.META_WAVELENGTH: self.data_df[cols.WAVELENGTH],
                cols.META_ABSORBANCE: meta_abs_mean,
            }
        )

    def irr1_sensitivity(self):
        """Return the PSS at ``irr1`` for every value in ``x_array``."""
        return np.array([
            self.irr1_pss(X)
            for X in self.x_array
        ])

    def irr2_sensitivity(self):
        """Return the PSS at ``irr2`` for every value in ``x_array``."""
        return np.array([
            self.irr2_pss(X)
            for X in self.x_array
        ])
    
    def fischer_result(self, X=1.0):
        """Return pair identifiers and calculated PSS values at one QY ratio."""
        irr1_pss, irr2_pss = self.pss_values(X)

        return {
            'irr1': self.irr1,
            'irr2': self.irr2,
            'X value': X,
            'irr1_pss_X': irr1_pss,
            'irr2_pss_X': irr2_pss,
        }

    def sensitivity_result(self):
        """Summarise single-pair sensitivity results.

        Returns
        -------
        dict
            Pair labels, quantum-yield-ratio values, PSS sensitivity arrays,
            standard Fischer results at ``X = 1``, and minimum and maximum finite
            PSS values for each irradiation wavelength.
        """
        irr1_pss = self.irr1_sensitivity()
        irr2_pss = self.irr2_sensitivity()

        return {
            'irr1': self.irr1,
            'irr2': self.irr2,
            'x_array': self.x_array,
            'irr1_pss': irr1_pss,
            'irr2_pss': irr2_pss,
            'irr1_pss_X1': self.irr1_pss(1.0),
            'irr2_pss_X1': self.irr2_pss(1.0),
            'irr1_min': np.nanmin(irr1_pss),
            'irr1_max': np.nanmax(irr1_pss),
            'irr2_min': np.nanmin(irr2_pss),
            'irr2_max': np.nanmax(irr2_pss),
        }


@dataclass
class FischerMultiPairAnalysis:
    """Perform multipair Fischer analysis across all irradiation pairs.

    One :class:`FischerPairAnalysis` object is constructed for every unique
    combination of irradiation wavelengths. These objects are reused when
    constructing PSS and metastate matrices at different quantum-yield ratios.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Validated and normalised UV–Vis data.
    irr_wls : list[str]
        Irradiation wavelength labels.
    darkmax_range : tuple[float, float]
        Wavelength range, in nm, used to locate the dark-spectrum absorbance
        maximum.
    x_array : numpy.ndarray
        Quantum-yield-ratio values used for sensitivity analysis.
    X : float, optional
        Default quantum-yield ratio used when matrix-building methods receive
        no explicit value.

    Attributes
    ----------
    pairs : list[FischerPairAnalysis]
        Analyses for every unique irradiation-wavelength pair.
    """
    data_df: pd.DataFrame
    irr_wls: list[str]
    darkmax_range: tuple[float, float]
    x_array: np.ndarray
    X: float = 1.0

    pairs: list[FischerPairAnalysis] = field(init=False)

    def __post_init__(self):
        """Construct an analysis object for every unique irradiation pair."""
        self.pairs = [
            FischerPairAnalysis(
                self.data_df,
                irr1,
                irr2,
                self.darkmax_range,
                self.x_array,
            )
            for irr1, irr2 in combinations(self.irr_wls, 2)
        ]

    def build_pss_matrix(self, X=None):
        """Construct the complete multipair PSS matrix at one ratio.

        For every Fischer pair, direct PSS values are calculated at the two pair
        wavelengths and PSS values are imputed for all remaining observed
        irradiation wavelengths.

        Parameters
        ----------
        X : float or None, optional
            Quantum-yield ratio. If ``None``, the instance default ``X`` is used.

        Returns
        -------
        pandas.DataFrame
            PSS values for every Fischer-pair and observed-irradiation combination.
        """
        if X is None:
            X = self.X

        rows = []

        for pair in self.pairs:
            result = pair.fischer_result(X)

            for observed_irr in self.irr_wls:
                if observed_irr == pair.irr1:
                    pss = result['irr1_pss_X']

                elif observed_irr == pair.irr2:
                    pss = result['irr2_pss_X']

                else:
                    pss = pair.impute_irr3_pss(
                        observed_irr,
                        X,
                    )

                rows.append(
                    {
                        cols.PAIR_IRR1: pair.irr1,
                        cols.PAIR_IRR2: pair.irr2,
                        cols.OBSERVED_IRR: observed_irr,
                        cols.PSS: pss,
                    }
                )

        return pd.DataFrame(rows)

    def build_meta_matrix(self, X=None):
        """Construct metastate spectra for every Fischer pair at one ratio.

        Parameters
        ----------
        X : float or None, optional
            Quantum-yield ratio. If ``None``, the instance default ``X`` is used.

        Returns
        -------
        pandas.DataFrame
            Long-form extrapolated metastate spectra for every Fischer pair.
        """
        if X is None:
            X = self.X

        meta_dfs = [
            pair.extrapolate_meta_spectrum(X)
            for pair in self.pairs
        ]

        meta_matrix = pd.concat(meta_dfs, ignore_index=True)

        return meta_matrix
    
    def build_matrices(self, X=None):
        """Construct the PSS and metastate matrices at one ratio.

        Parameters
        ----------
        X : float or None, optional
            Quantum-yield ratio. If ``None``, the instance default ``X`` is used.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Multipair PSS matrix and long-form metastate spectra matrix.
        """
        return (
            self.build_pss_matrix(X),
            self.build_meta_matrix(X),
        )

    def build_error_matrix(
        self,
        pss_matrix,
        final_pss_summary,
    ):
        """Build the directed direct-pair absolute-error matrix.

        Direct-pair PSS values at ``X = 1`` are compared with the final accepted
        multipair mean for the same target wavelength. Each row represents the
        target wavelength and each column represents its partner. Error is reported
        in percentage points.

        Parameters
        ----------
        pss_matrix : pandas.DataFrame
            Complete PSS matrix calculated at ``X = 1``.
        final_pss_summary : pandas.DataFrame
            Final accepted PSS summary for every observed irradiation wavelength.

        Returns
        -------
        tuple[list[str], numpy.ndarray]
            Ordered irradiation wavelengths and the square absolute-error matrix.
        """

        error_matrix = np.full(
            (len(self.irr_wls), len(self.irr_wls)),
            np.nan,
            dtype=float,
        )

        if pss_matrix.empty or final_pss_summary.empty:
            return list(self.irr_wls), error_matrix

        irr_index = {
            irradiation: index
            for index, irradiation in enumerate(self.irr_wls)
        }

        final_mean_lookup = (
            final_pss_summary
            .set_index(cols.OBSERVED_IRR)[cols.MEAN_PSS]
        )

        direct_rows = pss_matrix.loc[
            pss_matrix[cols.OBSERVED_IRR].eq(
                pss_matrix[cols.PAIR_IRR1]
            )
            | pss_matrix[cols.OBSERVED_IRR].eq(
                pss_matrix[cols.PAIR_IRR2]
            )
        ]

        for _, row in direct_rows.iterrows():
            irr1 = row[cols.PAIR_IRR1]
            irr2 = row[cols.PAIR_IRR2]
            target = row[cols.OBSERVED_IRR]
            pair_pss = float(row[cols.PSS])

            if target not in final_mean_lookup.index:
                continue

            final_mean = float(
                final_mean_lookup.loc[target]
            )

            if not (
                np.isfinite(pair_pss)
                and np.isfinite(final_mean)
            ):
                continue

            partner = (
                irr2
                if target == irr1
                else irr1
            )

            error_matrix[
                irr_index[target],
                irr_index[partner],
            ] = 100.0 * abs(
                pair_pss - final_mean
            )

        return list(self.irr_wls), error_matrix


    def build_sensitivity_matrix(self):
        """Build the directed direct-pair QY-sensitivity matrix.

        Each row represents the target irradiation wavelength and each column
        represents its partner wavelength. Sensitivity is the range of finite PSS
        values calculated across ``x_array`` and is reported in percentage points.

        Returns
        -------
        tuple[list[str], numpy.ndarray]
            Ordered irradiation wavelengths and the square sensitivity matrix.
        """

        irr_index = {
            irradiation: index
            for index, irradiation in enumerate(self.irr_wls)
        }

        sensitivity_matrix = np.full(
            (len(self.irr_wls), len(self.irr_wls)),
            np.nan,
            dtype=float,
        )

        for pair in self.pairs:
            for target, partner, sensitivity_values in (
                (
                    pair.irr1,
                    pair.irr2,
                    pair.irr1_sensitivity(),
                ),
                (
                    pair.irr2,
                    pair.irr1,
                    pair.irr2_sensitivity(),
                ),
            ):
                sensitivity_values = np.asarray(
                    sensitivity_values,
                    dtype=float,
                )

                finite_values = sensitivity_values[
                    np.isfinite(sensitivity_values)
                ]

                if finite_values.size == 0:
                    continue

                sensitivity = (
                    np.max(finite_values)
                    - np.min(finite_values)
                )

                sensitivity_matrix[
                    irr_index[target],
                    irr_index[partner],
                ] = 100.0 * float(sensitivity)

        return list(self.irr_wls), sensitivity_matrix


def qyratio_array(qyratio_range, n_points = N_QY_POINTS):
    """Generate evenly spaced quantum-yield-ratio values.

    Parameters
    ----------
    qyratio_range : iterable of float
        Lower and upper quantum-yield-ratio limits.
    n_points : int, optional
        Number of sampled values.

    Returns
    -------
    numpy.ndarray
        Evenly spaced quantum-yield-ratio values including both limits.
    """

    lower, upper = [float(value) for value in qyratio_range]
    return np.linspace(lower, upper, int(n_points))


def extract_abs_data(df, irr1, irr2, idx_darkmax):
    """Extract absorbance values required for one Fischer calculation.

    Experimental spectral rows nearest to ``irr1`` and ``irr2`` are used for
    the irradiation-wavelength absorbances.

    Parameters
    ----------
    df : pandas.DataFrame
        Validated and normalised UV–Vis data.
    irr1, irr2 : str
        Irradiation wavelength labels defining the Fischer pair.
    idx_darkmax : int
        DataFrame index of the selected dark-spectrum absorbance maximum.

    Returns
    -------
    AbsorbanceData
        Absorbance values required by the Fischer equations.
    """
    wl_data = df[cols.WAVELENGTH]
    
    irr1_index = utils.nearest_index(wl_data, float(irr1))
    irr2_index = utils.nearest_index(wl_data, float(irr2))

    return AbsorbanceData(
        dark_irr1=float(df[cols.DARK].iat[irr1_index]),
        dark_irr2=float(df[cols.DARK].iat[irr2_index]),
        dark_max=float(df.at[idx_darkmax, cols.DARK]),
        irr1_max=float(df.at[idx_darkmax, irr1]),
        irr2_max=float(df.at[idx_darkmax, irr2]),
        irr1_at_irr1=float(df[irr1].iat[irr1_index]),
        irr2_at_irr2=float(df[irr2].iat[irr2_index]),
    )


def best_irr(df, irr_wls, darkmax_range):
    """Select the irradiation spectrum with the largest dark-maximum change.

    Parameters
    ----------
    df : pandas.DataFrame
        Validated and normalised UV–Vis data.
    irr_wls : list[str]
        Irradiation wavelength labels.
    darkmax_range : tuple[float, float]
        Wavelength range, in nm, used to locate the dark-spectrum absorbance
        maximum.

    Returns
    -------
    str
        Irradiation wavelength label with the largest absolute absorbance
        difference from the dark spectrum at its selected maximum.
    """

    idx_darkmax, _ = utils.dark_lambda_max(df, darkmax_range)

    irr_diffs = {
        irr_wl: abs(
            float(df.at[idx_darkmax, cols.DARK]) - float(df.at[idx_darkmax, irr_wl])
        )
        for irr_wl in irr_wls
    }

    bestirr = max(irr_diffs, key=irr_diffs.get)

    return str(bestirr)


def multifischer(abs_data, X, target='irr1'):
    """Calculate the PSS for one wavelength of a Fischer pair.

    The standard Fischer expression is used when the effective
    quantum-yield ratio equals one. Otherwise, the continuous solution of the
    MultiFischer quadratic is selected. For ``target='irr2'``, the pair ordering is
    reversed and the reciprocal quantum-yield ratio is used.

    Parameters
    ----------
    abs_data : AbsorbanceData
        Absorbance values for the Fischer pair in ``irr1``, ``irr2`` order.
    X : float
        Quantum-yield ratio for ``irr1`` relative to ``irr2``.
    target : {'irr1', 'irr2'}, optional
        Pair wavelength for which the PSS is calculated.

    Returns
    -------
    float
        Calculated PSS, or ``numpy.nan`` when the calculation is undefined or
        the quadratic discriminant is negative.

    Raises
    ------
    ValueError
        If ``target`` is neither ``"irr1"`` nor ``"irr2"``.
    """
    if target == 'irr1':
        calc_abs_data = abs_data
        calc_X = X

    elif target == 'irr2':
        calc_abs_data = abs_data.reversed()
        calc_X = 1.0 / X

    else:
        raise ValueError('target must be "irr1" or "irr2".')
    
    a_dark_irr1 = calc_abs_data.dark_irr1
    a_dark_irr2 = calc_abs_data.dark_irr2
    d1 = calc_abs_data.d1
    d2 = calc_abs_data.d2

    # The Fischer equations normalise absorbance changes by the dark
    # absorbance at both irradiation wavelengths. Exact zeros therefore make
    # this pair undefined, although zero absorbance is still valid input data.
    if a_dark_irr1 == 0.0 or a_dark_irr2 == 0.0:
        return np.nan

    denom_n = calc_abs_data.irr1_max - calc_abs_data.dark_max

    if denom_n == 0:
        return np.nan

    n_value = (calc_abs_data.irr2_max - calc_abs_data.dark_max) / denom_n

    a_coef = (calc_X * n_value) - n_value
    b_coef = 1.0 - (n_value * d1 / a_dark_irr1) - (calc_X * n_value) + (calc_X * d2 / a_dark_irr2)
    c_coef = (d1 / a_dark_irr1) - (calc_X * d2 / a_dark_irr2)
    b_1_coef = 1.0 - (n_value * d1 / a_dark_irr1) - n_value + (d2 / a_dark_irr2) # b_coef(X=1)

    if calc_X == 1: # Fischer PSS: MultiFischer PSS equation would return 0/0 = undefined at X=1 so use standard Fischer
        numerator = (d2 / a_dark_irr2) - (d1 / a_dark_irr1)
        denominator = 1.0 + (d2 / a_dark_irr2) - n_value * (1.0 + (d1 / a_dark_irr1))
        
        if denominator == 0:
            return np.nan

        pss = numerator / denominator

        return pss

    pss = solve_quadratic_continuous_root(a_coef, b_coef, c_coef, b_1_coef)

    return pss


def solve_quadratic_continuous_root(a, b, c, b_1):
    """Select the MultiFischer quadratic root that is continuous at ``X = 1``.

    At ``X = 1``, the quadratic coefficient tends to zero. The sign of the
    linear coefficient at this limit, ``b_1``, determines which quadratic root
    has a numerator that also tends to zero and therefore remains continuous.

    Parameters
    ----------
    a, b, c : float
        Quadratic coefficients.
    b_1 : float
        Value of the linear coefficient at ``X = 1``.

    Returns
    -------
    float
        Continuous quadratic root, the direct linear solution if the quadratic
        coefficient is zero, or ``numpy.nan`` if no valid continuous solution
        is available.
    """

    # If the standard Fischer denominator is zero, the root that should be
    # continuous at X = 1 is undefined and neither quadratic branch is valid.
    if b_1 == 0.0:
        return np.nan

    # Some otherwise valid spectra reduce the quadratic to a linear equation.
    # Solve that equation directly instead of dividing by 2a.
    if a == 0.0:
        if b == 0.0:
            return np.nan
        return -c / b

    discr = b**2 - 4.0 * a * c
    if discr < 0:
        return np.nan

    if b_1 > 0:
        return (-b + np.sqrt(discr)) / (2.0 * a)
    else:
        return (-b - np.sqrt(discr)) / (2.0 * a)


def summarise_accepted_pss(
    accepted_pss_matrix,
    observed_irrs=None,
):
    """Summarise PSS values passing every applied filter.

    The reported mean and sample standard deviation are calculated from all
    fully accepted values, not from the HDI core used to define the HDI
    acceptance bounds.

    Parameters
    ----------
    accepted_pss_matrix : pandas.DataFrame
        PSS matrix containing only rows that pass every applied filter.
    observed_irrs : iterable of str or None, optional
        Complete ordered set of observed irradiation wavelengths. If supplied,
        wavelengths with no accepted values are retained with a count of zero
        and missing mean and standard deviation.

    Returns
    -------
    pandas.DataFrame
        Mean PSS, sample standard deviation, and accepted-value count for each
        observed irradiation wavelength.
    """

    final_summary = (
        accepted_pss_matrix
        .groupby(cols.OBSERVED_IRR, sort=False)[cols.PSS]
        .agg(['mean', 'std', 'size'])
        .rename(
            columns={
                'mean': 'Mean PSS',
                'std': 'Standard deviation',
                'size': 'Number of accepted values',
            }
        )
    )

    if observed_irrs is not None:
        final_summary = final_summary.reindex(
            pd.Index(
                observed_irrs,
                name=cols.OBSERVED_IRR,
            )
        )

    final_summary['Number of accepted values'] = (
        final_summary['Number of accepted values']
        .fillna(0)
        .astype(int)
    )

    return final_summary.reset_index()
    

def highest_pss_metastate_spectrum(
    data_df,
    final_pss_summary,
):
    """Extrapolate the metastate spectrum from the highest-PSS irradiation.

    The irradiation wavelength with the highest final mean PSS is selected.
    The central metastate spectrum is calculated using its mean PSS, while the
    lower and upper absorbance bounds are obtained by repeating the
    extrapolation at mean PSS minus and plus one sample standard deviation.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Validated and normalised UV-Vis data.
    final_pss_summary : pandas.DataFrame
        Final mean PSS and sample standard deviation for each irradiation
        wavelength.

    Returns
    -------
    tuple[str, pandas.DataFrame]
        Selected irradiation wavelength and extrapolated metastate spectrum
        containing the central, lower, and upper absorbance values.

    Raises
    -------
    ValueError
        If the selected mean PSS or either one standard deviation PSS bound is
        zero and therefore cannot be used as an extrapolation denominator.
    """

    highest_pss_index = final_pss_summary[
        cols.MEAN_PSS
    ].idxmax()

    highest_pss_row = final_pss_summary.loc[
        highest_pss_index
    ]

    irradiation = highest_pss_row[
        cols.OBSERVED_IRR
    ]

    mean_pss = float(
        highest_pss_row[cols.MEAN_PSS]
    )

    std_pss = float(
        highest_pss_row[cols.PSS_STD]
    )

    pss_denominators = {
        'mean PSS': mean_pss,
        'mean PSS minus one standard deviation': mean_pss - std_pss,
        'mean PSS plus one standard deviation': mean_pss + std_pss,
    }

    zero_denominators = [
        label
        for label, value in pss_denominators.items()
        if value == 0.0
    ]

    if zero_denominators:
        raise ValueError(
            'Cannot extrapolate the metastate spectrum from '
            f'{irradiation} nm because {", ".join(zero_denominators)} is zero.'
        )

    dark_absorbance = data_df[
        cols.DARK
    ].to_numpy(dtype=float)

    absorbance_change = (
        data_df[irradiation].to_numpy(dtype=float)
        - dark_absorbance
    )

    mean_metastate = (
        dark_absorbance
        + absorbance_change / mean_pss
    )

    metastate_at_lower_pss = (
        dark_absorbance
        + absorbance_change
        / (mean_pss - std_pss)
    )

    metastate_at_upper_pss = (
        dark_absorbance
        + absorbance_change
        / (mean_pss + std_pss)
    )

    # Which PSS limit produces the lower absorbance can change across the
    # spectrum, so order the two curves independently at every wavelength
    lower_metastate = np.minimum(
        metastate_at_lower_pss,
        metastate_at_upper_pss,
    )

    upper_metastate = np.maximum(
        metastate_at_lower_pss,
        metastate_at_upper_pss,
    )

    metastate_spectrum = pd.DataFrame(
        {
            cols.META_WAVELENGTH: data_df[
                cols.WAVELENGTH
            ].to_numpy(dtype=float),
            cols.MEAN_META_ABSORBANCE: mean_metastate,
            cols.LOWER_META_ABSORBANCE: lower_metastate,
            cols.UPPER_META_ABSORBANCE: upper_metastate,
        }
    )

    return irradiation, metastate_spectrum


def singlepair(
    data_df,
    irr1,
    irr2,
    darkmax_range,
    qyratio_range,
):
    """Run quantum-yield-ratio sensitivity analysis for one Fischer pair.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Validated and normalised UV–Vis data.
    irr1, irr2 : str
        Irradiation wavelength labels defining the Fischer pair.
    darkmax_range : tuple[float, float]
        Wavelength range, in nm, used to locate the dark-spectrum absorbance
        maximum.
    qyratio_range : tuple[float, float]
        Lower and upper quantum-yield-ratio limits.

    Returns
    -------
    dict
        Single-pair sensitivity results, including the sampled ratios, PSS
        arrays, standard Fischer values, and sensitivity extrema.
    """
    outlog.log(f'  Quantum yield ratio (X) range: {qyratio_range[0]} - {qyratio_range[1]}')
    x_array = qyratio_array(qyratio_range)
    outlog.logln(f'  Generated X array with {N_QY_POINTS} datapoints.')

    pair = FischerPairAnalysis(
        data_df=data_df,
        irr1=irr1,
        irr2=irr2,
        darkmax_range=darkmax_range,
        x_array=x_array,
    )

    results = pair.sensitivity_result()

    return results


def multipair(
    data_df,
    irr_wls,
    darkmax_range,
    qyratio_range,
):
    """Run standard multipair Fischer analysis at ``X = 1``.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Validated and normalised UV-Vis data.
    irr_wls : list[str]
        Irradiation wavelength labels.
    darkmax_range : tuple[float, float]
        Wavelength range, in nm, used to locate the dark-spectrum absorbance
        maximum.
    qyratio_range : tuple[float, float]
        Quantum-yield-ratio limits used to initialise the associated pair
        analyses.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Complete PSS matrix and long-form metastate spectra matrix calculated at
        ``X = 1``.
    """
    x_array = qyratio_array(qyratio_range)

    analysis = FischerMultiPairAnalysis(
        data_df=data_df,
        irr_wls=irr_wls,
        darkmax_range=darkmax_range,
        x_array=x_array,
        X=1.0,
    )

    pss_matrix, meta_matrix = analysis.build_matrices()

    return pss_matrix, meta_matrix
