"""Multipair PSS filtering tools."""

from dataclasses import dataclass, field
import numpy as np
import pandas as pd

from . import columns as cols

@dataclass
class PSSFilterPipeline:
    """Apply physical and statistical filters to a multipair PSS matrix.

    The pipeline retains all input rows and adds Boolean columns describing the
    result of each applied filter. Filters are applied in sequence, allowing
    later filters to use the results of earlier filters.

    Parameters
    ----------
    pss_matrix : pandas.DataFrame
        Multipair PSS matrix containing Fischer-pair identifiers, observed
        irradiation wavelengths, and calculated PSS values.
    meta_matrix : pandas.DataFrame
        Long-form extrapolated metastate spectra for each Fischer pair.
    pss_range : tuple[float, float]
        Inclusive range of accepted PSS values, expressed as fractions.
    meta_range : tuple[float, float]
        Inclusive wavelength range, in nm, used for the metastate filter.
    meta_min : float
        Minimum permitted extrapolated metastate absorbance within
        ``meta_range``.
    hdi : float
        Fraction of candidate PSS values included in the empirical HDI core.
    hdi_sigma : float
        Number of HDI-core sample standard deviations used to extend the HDI
        acceptance bounds.

    Attributes
    ----------
    filters : list[str]
        Names of the filter-result columns that have been applied.
    hdi_diagnostics : pandas.DataFrame
        HDI statistics for each observed irradiation wavelength.
    """

    pss_matrix: pd.DataFrame
    meta_matrix: pd.DataFrame
    pss_range: tuple[float, float]
    meta_range: tuple[float, float]
    meta_min: float
    hdi: float
    hdi_sigma: float

    filters: list[str] = field(default_factory=list, init=False)

    # Store one HDI diagnostic row for each observed irradiation wavelength
    hdi_diagnostics: pd.DataFrame = field(
        default_factory=pd.DataFrame,
        init=False,
    )

    def __post_init__(self):
        """Copy the input matrices so that filtering does not modify the originals."""

        self.pss_matrix = self.pss_matrix.copy()
        self.meta_matrix = self.meta_matrix.copy()

    def add_filter(self, filter_name):
        """Record an applied filter if it has not already been recorded."""

        if filter_name not in self.filters:
            self.filters.append(filter_name)


    def apply_pss_range_filter(self):
        """Apply the pair-level PSS-range filter.

        Each individual PSS value is tested against ``pss_range``. A Fischer pair
        passes only if its PSS values for every observed irradiation wavelength are
        finite and within the inclusive range. The pair-level result is assigned to
        every row belonging to that pair.

        Returns
        -------
        PSSFilterPipeline
            The current pipeline, allowing chained method calls.
        """

        lower, upper = self.pss_range

        row_pass = (
            np.isfinite(self.pss_matrix[cols.PSS])
            & self.pss_matrix[cols.PSS].between(
                lower,
                upper,
                inclusive='both',
            )
        )

        # Require every observed-irradiation row in a Fischer pair to pass
        pair_pass = row_pass.groupby(
            [
                self.pss_matrix[cols.PAIR_IRR1],
                self.pss_matrix[cols.PAIR_IRR2],
            ]
        ).transform('all')

        self.pss_matrix[cols.PSS_FILTER] = pair_pass.astype(bool)
        self.add_filter(cols.PSS_FILTER)

        return self


    def apply_metastate_filter(self):
        """Apply the pair-level metastate-absorbance filter.

        For each Fischer pair, the minimum extrapolated metastate absorbance within
        ``meta_range`` is calculated. The pair passes if this minimum is greater
        than or equal to ``meta_min``. The pair-level result is assigned to every
        PSS row belonging to that pair.

        Returns
        -------
        PSSFilterPipeline
            The current pipeline, allowing chained method calls.
        """

        meta_min_wl, meta_max_wl = self.meta_range

        meta_in_wl_range = self.meta_matrix[
            self.meta_matrix[cols.META_WAVELENGTH].between(
                meta_min_wl,
                meta_max_wl,
                inclusive='both',
            )
        ]

        meta_filter_summary = (
            meta_in_wl_range
            .groupby(
                [cols.PAIR_IRR1, cols.PAIR_IRR2],
                as_index=False,
            )[cols.META_ABSORBANCE]
            .min()
            .rename(
                columns={
                    cols.META_ABSORBANCE: cols.MIN_META_ABSORBANCE,
                }
            )
        )

        meta_filter_summary[cols.METASTATE_FILTER] = (
            meta_filter_summary[cols.MIN_META_ABSORBANCE]
            >= self.meta_min
        )

        self.pss_matrix = self.pss_matrix.drop(
            columns=[
                column for column in (cols.MIN_META_ABSORBANCE, cols.METASTATE_FILTER)
                if column in self.pss_matrix.columns
            ]
        )

        self.pss_matrix = self.pss_matrix.merge(
            meta_filter_summary,
            on=[cols.PAIR_IRR1, cols.PAIR_IRR2],
            how='left',
        )

        self.pss_matrix[cols.METASTATE_FILTER] = (
            self.pss_matrix[cols.METASTATE_FILTER]
            .fillna(False)
            .astype(bool)
        )

        self.add_filter(cols.METASTATE_FILTER)

        return self


    def apply_hdi_filter(self):
        """Apply an HDI filter independently to each observed irradiation wavelength.

        Candidate PSS values are selected separately for each observed irradiation
        wavelength from rows passing every previously applied filter. The narrowest
        contiguous interval containing at least the requested ``hdi`` fraction is
        used as the empirical HDI core.

        Extended acceptance bounds are defined from the HDI-core mean and sample
        standard deviation as:

            lower = core_mean - hdi_sigma * core_std
            upper = core_mean + hdi_sigma * core_std

        Candidate values within the inclusive bounds pass the HDI filter. The
        HDI-core statistics are used only to define these bounds; final reported PSS
        statistics are calculated separately from all fully accepted rows.

        If no candidates are available, diagnostic values are recorded as missing.
        If the HDI core contains one value, its sample standard deviation and
        acceptance bounds are missing, so no rows pass the HDI filter.

        Returns
        -------
        PSSFilterPipeline
            The current pipeline, allowing chained method calls.
        """

        # Restrict HDI candidates to rows passing every earlier filter
        previous_filters = [
            filter_name
            for filter_name in self.filters
            if filter_name != cols.HDI_FILTER
        ]

        if previous_filters:
            previous_filter_mask = (
                self.pss_matrix[previous_filters]
                .all(axis=1)
            )
        else:
            # Allow HDI filtering when no earlier filters have been applied
            previous_filter_mask = pd.Series(
                True,
                index=self.pss_matrix.index,
                dtype=bool,
            )
        
        # Remove NaN or infinite PSS values
        finite_pss_mask = pd.Series(
            np.isfinite(
                self.pss_matrix[cols.PSS].to_numpy(dtype=float)
            ),
            index=self.pss_matrix.index,
            dtype=bool,
        )

        # Begin with every row rejected. A row is accepted only if it passes
        # the HDI bounds for its own observed irradiation wavelength
        self.pss_matrix[cols.HDI_FILTER] = False

        diagnostic_rows = []

        # Give each observed irradiation wavelength its own HDI core and
        # extended acceptance bounds
        for observed_irr in self.pss_matrix[cols.OBSERVED_IRR].drop_duplicates():

            irradiation_mask = self.pss_matrix[cols.OBSERVED_IRR].eq(observed_irr)

            candidate_mask = (
                irradiation_mask
                & previous_filter_mask
                & finite_pss_mask
            )

            candidate_values = np.sort(
                self.pss_matrix.loc[
                    candidate_mask,
                    cols.PSS,
                ].to_numpy(dtype=float)
            )

            n_candidates = len(candidate_values)

            # Retain a diagnostic row even when no values survive the earlier
            # filters. This ensures every observed irradiation wavelength can
            # still appear in the exported csv
            if n_candidates == 0:
                diagnostic_rows.append(
                    {
                        cols.OBSERVED_IRR: observed_irr,
                        'Number of candidates': 0,
                        'HDI fraction': self.hdi,
                        'HDI core size': 0,
                        'HDI core minimum': np.nan,
                        'HDI core maximum': np.nan,
                        'HDI core mean': np.nan,
                        'HDI core standard deviation': np.nan,
                        'HDI sigma': self.hdi_sigma,
                        'HDI lower bound': np.nan,
                        'HDI upper bound': np.nan,
                        'HDI pass': 0,
                        'HDI fail': 0,
                    }
                )
                continue

            # ceil ensures the HDI core contains at least the requested fraction
            hdi_n = max(1, int(np.ceil(self.hdi * n_candidates)))

            # Because candidate_values is sorted, each possible contiguous HDI
            # core is a window of hdi_n values. Its width is the final value
            # minus the first value
            window_widths = (
                candidate_values[hdi_n - 1:]
                - candidate_values[:n_candidates - hdi_n + 1]
            )

            # Select the narrowest window as the highest-density core
            hdi_start = int(np.argmin(window_widths))
            hdi_core = candidate_values[
                hdi_start:hdi_start + hdi_n
            ]

            hdi_core_mean = np.mean(hdi_core)

            # Use the sample standard deviation for multi-value cores. A
            # single-value core has an undefined (NaN) standard deviation
            hdi_core_std = (
                np.std(hdi_core, ddof=1)
                if hdi_n > 1
                else np.nan
            )

            # These filtering statistics are distinct from the final reported
            # mean and standard deviation
            hdi_lower = (
                hdi_core_mean
                - self.hdi_sigma * hdi_core_std
            )
            hdi_upper = (
                hdi_core_mean
                + self.hdi_sigma * hdi_core_std
            )

            # Only valid candidates can pass the HDI filter, so rows rejected
            # by an earlier filter remain False in "HDI filter"
            wavelength_pass_mask = (
                candidate_mask
                & self.pss_matrix[cols.PSS].between(
                    hdi_lower,
                    hdi_upper,
                    inclusive='both',
                )
            )

            self.pss_matrix.loc[
                wavelength_pass_mask,
                cols.HDI_FILTER,
            ] = True

            n_pass = int(wavelength_pass_mask.sum())

            diagnostic_rows.append(
                {
                    cols.OBSERVED_IRR: observed_irr,
                    'Number of candidates': int(n_candidates),
                    'HDI fraction': float(self.hdi),
                    'HDI core size': int(hdi_n),
                    'HDI core minimum': float(hdi_core[0]),
                    'HDI core maximum': float(hdi_core[-1]),
                    'HDI core mean': hdi_core_mean,
                    'HDI core standard deviation': hdi_core_std,
                    'HDI sigma': float(self.hdi_sigma),
                    'HDI lower bound': hdi_lower,
                    'HDI upper bound': hdi_upper,
                    'HDI pass': n_pass,
                    'HDI fail': int(n_candidates - n_pass),
                }
            )

        self.hdi_diagnostics = pd.DataFrame(diagnostic_rows)

        self.add_filter(cols.HDI_FILTER)

        return self


    def all_filter_mask(self):
        """Return a Boolean mask identifying rows that pass every applied filter."""

        if not self.filters:
            return pd.Series(True, index=self.pss_matrix.index)

        return self.pss_matrix[self.filters].all(axis=1)

    def accepted_matrix(self):
        """Return a copy of the rows that pass every applied filter."""

        return self.pss_matrix.loc[self.all_filter_mask()].copy()

    def filter_counts(self):
        """Return pass and fail counts for individual and combined filters."""

        counts = {
            'Total rows': int(len(self.pss_matrix)),
        }

        for column in self.filters:
            counts[f'{column} pass'] = int(self.pss_matrix[column].sum())
            counts[f'{column} fail'] = int((~self.pss_matrix[column]).sum())

        counts['All filters pass'] = int(self.all_filter_mask().sum())
        counts['All filters fail'] = int((~self.all_filter_mask()).sum())

        return counts