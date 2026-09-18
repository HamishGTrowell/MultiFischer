"""Plotting tools for MultiFischer results."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm

from . import outlog
from . import columns as cols
from . import utils

plt.style.use('seaborn-v0_8-colorblind')

GLOBAL_PLOT_STYLE_UPDATES = {
#    'font.size': 18,
#    'axes.labelsize': 18,
#    'axes.titlesize': 20,
#    'legend.fontsize': 16,
}

plt.rcParams.update(GLOBAL_PLOT_STYLE_UPDATES)


def irradiation_color_map(irradiations):
    """Return a consistent colour for each irradiation wavelength."""

    irradiations = sorted(
        set(irradiations),
        key=utils.wavelength_value_to_float,
    )

    if not irradiations:
        return {}

    cmap = plt.get_cmap(
        'viridis',
        len(irradiations),
    )

    return {
        irradiation: cmap(index)
        for index, irradiation in enumerate(irradiations)
    }


def violin(
    pss_matrix,
    observed_irrs,
    metastate,
    pss_range,
    out,
):
    """Plot PSS distributions and filter outcomes for selected wavelengths.

    The violin shapes show the complete finite PSS distribution for each
    selected observed irradiation wavelength. Individual points are coloured
    according to the first filtering stage at which they were rejected, or
    shown as accepted if they passed every filter.

    Parameters
    ----------
    pss_matrix : pandas.DataFrame
        Complete PSS matrix containing the PSS values and Boolean filter-result
        columns.
    observed_irrs : iterable of str
        Observed irradiation wavelengths to include in the plot.
    metastate : str
        Label of the metastable isomer used in the y-axis label.
    pss_range : tuple[float, float]
        Accepted PSS range, used to determine the displayed y-axis limits.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path or None
        Path to the saved PNG file, or ``None`` if no finite PSS values are
        available for the selected wavelengths.
    """

    observed_irrs = list(observed_irrs)

    plot_data = pss_matrix.loc[
        pss_matrix[cols.OBSERVED_IRR].isin(observed_irrs)
        & np.isfinite(pss_matrix[cols.PSS])
        & pss_matrix[cols.PSS].between(
            pss_range[0] - 0.1,
            pss_range[1] + 0.1,
            inclusive='both',
        )
    ].copy()

    if plot_data.empty:
        outlog.warning(
            'No finite PSS values available for the requested violin plot.'
        )
        return None

    outcome_col = 'Filter outcome'

    plot_data[outcome_col] = np.select(
        [
            ~plot_data[cols.PSS_FILTER],
            (
                plot_data[cols.PSS_FILTER]
                & ~plot_data[cols.METASTATE_FILTER]
            ),
            (
                plot_data[cols.PSS_FILTER]
                & plot_data[cols.METASTATE_FILTER]
                & ~plot_data[cols.HDI_FILTER]
            ),
        ],
        [
            'Failed PSS range',
            'Failed metastate',
            'Failed HDI',
        ],
        default='Accepted',
    )

    # Preserve the wavelength order requested by the user rather than allowing
    # Seaborn to sort the category labels alphabetically.
    plot_data[cols.OBSERVED_IRR] = pd.Categorical(
        plot_data[cols.OBSERVED_IRR],
        categories=observed_irrs,
        ordered=True,
    )

    outcome_order = [
        'Failed PSS range',
        'Failed metastate',
        'Failed HDI',
        'Accepted',
    ]

    palette = {
        'Failed PSS range': '#d73027',
        'Failed metastate': '#f39c12',
        'Failed HDI': '#756bb1',
        'Accepted': '#111111',
    }

    n_irrs = len(observed_irrs)
    figure_width = 4 if n_irrs == 1 else max(6, 1.2 * n_irrs)

    fig, ax = plt.subplots(
        figsize=(figure_width, 5),
    )

    sns.violinplot(
        data=plot_data,
        x=cols.OBSERVED_IRR,
        y=cols.PSS,
        order=observed_irrs,
        color='lightblue',
        inner=None,
        cut=2,
        width=0.8 if n_irrs == 1 else 0.95,
        ax=ax,
    )

    sns.swarmplot(
        data=plot_data,
        x=cols.OBSERVED_IRR,
        y=cols.PSS,
        hue=outcome_col,
        order=observed_irrs,
        hue_order=outcome_order,
        palette=palette,
        size=5,
        ax=ax,
    )

    ax.yaxis.set_major_formatter(
        PercentFormatter(xmax=1.0)
    )

    if n_irrs == 1:
        ax.set_xlabel('')
        ax.set_xticks(
            ax.get_xticks(),
            labels=[f'{observed_irrs[0]} nm'],
        )
    else:
        ax.set_xlabel('Irradiation wavelength (nm)')
        ax.set_xticks(
            ax.get_xticks(),
            labels=observed_irrs,
        )

    ax.set_ylabel(
        f'PSS distribution (% $\\mathit{{{metastate}}}$)'
    )

    ax.set_ylim(
        float(pss_range[0]) - 0.1,
        float(pss_range[1]) + 0.1,
    )

    ax.grid(
        True,
        axis='y',
        linestyle='--',
        alpha=0.4,
    )

    # Seaborn normally includes only outcomes present in the plotted data.
    # This avoids displaying unused categories while preserving their order.
    handles, labels = ax.get_legend_handles_labels()

    if handles:
        ax.legend(
            handles,
            labels,
            title='Filter outcome',
            bbox_to_anchor=(1.02, 1.0),
            loc='upper left',
            borderaxespad=0.0,
        )

    if n_irrs == 1:
        filename = f'pss_distribution_{observed_irrs[0]}nm.png'
    elif observed_irrs == list(
        pss_matrix[cols.OBSERVED_IRR].drop_duplicates()
    ):
        filename = 'pss_distribution_all_irrs.png'
    else:
        wavelength_label = '_'.join(observed_irrs)
        filename = f'pss_distribution_{wavelength_label}nm.png'

    out_path = Path(out) / filename

    fig.tight_layout()
    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)

    return out_path


def spread_label_positions(values, lower, upper, min_separation):
    """Adjust sorted label positions to reduce vertical overlap."""

    if not values:
        return []

    positions = [float(values[0])]

    for value in values[1:]:
        positions.append(
            max(
                float(value),
                positions[-1] + min_separation,
            )
        )

    # Move all labels down if the final label exceeds the upper bound.
    overflow = positions[-1] - upper

    if overflow > 0.0:
        positions = [
            position - overflow
            for position in positions
        ]

    # Work backwards to restore the minimum separation if moving the labels
    # down caused earlier labels to overlap.
    for index in range(len(positions) - 2, -1, -1):
        positions[index] = min(
            positions[index],
            positions[index + 1] - min_separation,
        )

    # Move all labels up if the first label now falls below the lower bound.
    underflow = lower - positions[0]

    if underflow > 0.0:
        positions = [
            position + underflow
            for position in positions
        ]

    return positions


def multipair_sensitivity(
    sensitivity_summary,
    metastate,
    out,
):
    """Plot filtered mean PSS against QY ratio for each irradiation wavelength.

    Each curve shows the final mean PSS calculated from values passing all
    filters at each quantum-yield ratio. Missing means are retained as gaps in
    the curves rather than being connected across unavailable regions.

    Parameters
    ----------
    sensitivity_summary : pandas.DataFrame
        Long-form multipair sensitivity results containing the quantum-yield
        ratio, observed irradiation wavelength, final mean PSS, and accepted
        value count.
    metastate : str
        Label of the metastable isomer used in the y-axis label.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path or None
        Path to the saved PNG file, or ``None`` if no finite sensitivity
        results are available.
    """

    out_path = Path(out) / 'multipair_qy_sensitivity.png'

    if sensitivity_summary.empty:
        outlog.warning(
            'Multipair QY-sensitivity plot skipped because the sensitivity '
            'summary is empty.'
        )
        return None

    required_columns = {
        cols.QY_RATIO,
        cols.OBSERVED_IRR,
        cols.MEAN_PSS,
        cols.PSS_STD,
    }

    missing_columns = required_columns.difference(
        sensitivity_summary.columns
    )

    if missing_columns:
        raise ValueError(
            'Cannot create the multipair QY-sensitivity plot because the '
            f'following columns are missing: {sorted(missing_columns)}.'
        )

    plot_data = sensitivity_summary.copy()

    plot_data[cols.QY_RATIO] = pd.to_numeric(
        plot_data[cols.QY_RATIO],
        errors='coerce',
    )
    plot_data[cols.MEAN_PSS] = pd.to_numeric(
        plot_data[cols.MEAN_PSS],
        errors='coerce',
    )

    plot_data[cols.PSS_STD] = pd.to_numeric(
        plot_data[cols.PSS_STD],
        errors='coerce',
    )

    irradiations = sorted(
        plot_data[cols.OBSERVED_IRR]
        .dropna()
        .unique(),
        key=utils.wavelength_value_to_float,
    )

    color_map = irradiation_color_map(irradiations)

    fig, ax = plt.subplots(
        figsize=(8.8, 5.5),
    )

    endpoints = []

    for irradiation in irradiations:
        irradiation_data = (
            plot_data.loc[
                plot_data[cols.OBSERVED_IRR].eq(irradiation)
            ]
            .sort_values(cols.QY_RATIO)
        )

        finite_x = np.isfinite(
            irradiation_data[cols.QY_RATIO].to_numpy(dtype=float)
        )

        # Remove invalid X values, but retain NaN mean values so Matplotlib
        # leaves gaps where no PSS values passed all filters.
        irradiation_data = irradiation_data.loc[finite_x]

        if irradiation_data.empty:
            continue

        x_values = irradiation_data[
            cols.QY_RATIO
        ].to_numpy(dtype=float)

        mean_values = irradiation_data[
            cols.MEAN_PSS
        ].to_numpy(dtype=float)

        std_values = irradiation_data[
            cols.PSS_STD
        ].to_numpy(dtype=float)

        finite_means = np.isfinite(mean_values)

        if not np.any(finite_means):
            continue

        color = color_map[irradiation]

        lower_values = mean_values - std_values
        upper_values = mean_values + std_values

        # Draw the standard-deviation band first so the mean curve remains visible.
        # NaN standard deviations automatically create gaps in the band, including
        # points for which only one accepted PSS value was available.
        ax.fill_between(
            x_values,
            lower_values,
            upper_values,
            where=(
                np.isfinite(x_values)
                & np.isfinite(lower_values)
                & np.isfinite(upper_values)
            ),
            color=color,
            alpha=0.18,
            linewidth=0,
        )

        ax.plot(
            x_values,
            mean_values,
            color=color,
            linewidth=1.8,
        )

        # Label the last finite mean, even if the final sampled X value has no
        # accepted PSS results.
        last_finite_index = np.flatnonzero(finite_means)[-1]

        endpoints.append(
            {
                'x': float(x_values[last_finite_index]),
                'y': float(mean_values[last_finite_index]),
                'label': f'{irradiation} nm',
                'color': color,
            }
        )

    if not endpoints:
        plt.close(fig)
        outlog.warning(
            'Multipair QY-sensitivity plot skipped because no finite mean PSS '
            'values are available.'
        )
        return None

    ax.yaxis.set_major_formatter(
        PercentFormatter(xmax=1.0)
    )

    ax.set_xlabel(
        'Quantum-yield ratio (X)'
    )
    ax.set_ylabel(
        f'Final mean PSS (% $\\mathit{{{metastate}}}$)'
    )

    ax.grid(
        True,
        linestyle='--',
        alpha=0.4,
    )

    # Add space to the right of the curves for direct wavelength labels.
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()

    x_range = x_max - x_min
    y_range = y_max - y_min

    label_x = x_max + 0.04 * x_range

    ax.set_xlim(
        x_min,
        x_max + 0.22 * x_range,
    )

    endpoints.sort(
        key=lambda endpoint: endpoint['y']
    )

    label_positions = spread_label_positions(
        values=[
            endpoint['y']
            for endpoint in endpoints
        ],
        lower=y_min + 0.02 * y_range,
        upper=y_max - 0.02 * y_range,
        min_separation=0.04 * y_range,
    )

    for endpoint, label_y in zip(
        endpoints,
        label_positions,
    ):
        ax.plot(
            [endpoint['x'], label_x],
            [endpoint['y'], label_y],
            color=endpoint['color'],
            linewidth=0.8,
            alpha=0.75,
        )

        ax.text(
            label_x,
            label_y,
            endpoint['label'],
            color=endpoint['color'],
            fontsize=10,
            horizontalalignment='left',
            verticalalignment='center',
        )

    fig.tight_layout()
    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)

    return out_path


def plot_log_heatmap(
    matrix,
    irradiation_wavelengths,
    colorbar_label,
    filename,
    out,
):
    """Plot a positive numeric matrix using logarithmic colour scaling.

    Parameters
    ----------
    matrix : numpy.ndarray
        Square matrix containing positive values and optional missing cells.
    irradiation_wavelengths : iterable of str
        Irradiation wavelength labels defining both matrix axes.
    colorbar_label : str
        Label displayed beside the colour bar.
    filename : str
        Output PNG filename.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path or None
        Path to the saved PNG file, or ``None`` if no positive finite values
        are available.
    """

    matrix = np.asarray(
        matrix,
        dtype=float,
    )

    positive_mask = (
        np.isfinite(matrix)
        & (matrix > 0.0)
    )

    if not np.any(positive_mask):
        outlog.warning(
            f'{filename} skipped because the matrix contains no positive '
            'finite values.'
        )
        return None

    positive_values = matrix[positive_mask]

    fig, ax = plt.subplots(
        figsize=(7.5, 6.5),
    )

    sns.heatmap(
        matrix,
        cmap='viridis',
        norm=LogNorm(
            vmin=float(np.min(positive_values)),
            vmax=float(np.max(positive_values)),
        ),
        square=True,
        linewidths=0.1,
        xticklabels=irradiation_wavelengths,
        yticklabels=irradiation_wavelengths,
        cbar_kws={
            'label': colorbar_label,
        },
        ax=ax,
    )

    ax.set_xlabel(
        'Partner irradiation wavelength (nm)'
    )
    ax.set_ylabel(
        'Target irradiation wavelength (nm)'
    )

    ax.tick_params(
        axis='x',
        labelrotation=90,
    )
    ax.tick_params(
        axis='y',
        labelrotation=0,
    )

    out_path = Path(out) / filename

    fig.tight_layout()
    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)

    return out_path


def error_sensitivity_scatter(
    sensitivity_matrix,
    error_matrix,
    metastate,
    out,
):
    """Plot direct-pair sensitivity against error on logarithmic axes.

    A linear regression is fitted in base-10 logarithmic space:

        log10(error) = slope * log10(sensitivity) + intercept

    The regression line, equation, and coefficient of determination are shown
    on the plot. Only corresponding matrix cells containing positive finite
    sensitivity and error values are included.

    Parameters
    ----------
    sensitivity_matrix : numpy.ndarray
        Directed direct-pair sensitivity matrix in PSS percentage points.
    error_matrix : numpy.ndarray
        Directed direct-pair absolute-error matrix in PSS percentage points.
    metastate : str
        Label of the metastable isomer used in the axis labels.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path or None
        Path to the saved PNG file, or ``None`` if no valid paired values are
        available.
    """

    sensitivity = np.asarray(
        sensitivity_matrix,
        dtype=float,
    )
    error = np.asarray(
        error_matrix,
        dtype=float,
    )

    if sensitivity.shape != error.shape:
        raise ValueError(
            'Sensitivity and error matrices must have the same shape.'
        )

    sensitivity = sensitivity.ravel()
    error = error.ravel()

    valid_mask = (
        np.isfinite(sensitivity)
        & np.isfinite(error)
        & (sensitivity > 0.0)
        & (error > 0.0)
    )

    sensitivity = sensitivity[valid_mask]
    error = error[valid_mask]

    if sensitivity.size == 0:
        outlog.warning(
            'Error-sensitivity scatter plot skipped because no corresponding '
            'positive finite sensitivity and error values are available.'
        )
        return None

    fig, ax = plt.subplots(
        figsize=(6.5, 5.5),
    )

    ax.scatter(
        sensitivity,
        error,
        s=45,
        alpha=0.8,
    )

    ax.set_xscale('log')
    ax.set_yscale('log')

    # At least two distinct sensitivity values are required to fit a linear
    # regression in log space.
    if (
        sensitivity.size >= 2
        and np.unique(sensitivity).size >= 2
    ):
        log_sensitivity = np.log10(sensitivity)
        log_error = np.log10(error)

        slope, intercept = np.polyfit(
            log_sensitivity,
            log_error,
            deg=1,
        )

        predicted_log_error = (
            slope * log_sensitivity
            + intercept
        )

        residual_sum_squares = np.sum(
            (log_error - predicted_log_error) ** 2
        )
        total_sum_squares = np.sum(
            (log_error - np.mean(log_error)) ** 2
        )

        r_squared = (
            1.0
            - residual_sum_squares / total_sum_squares
            if total_sum_squares > 0.0
            else np.nan
        )

        trend_sensitivity = np.logspace(
            np.log10(np.min(sensitivity)),
            np.log10(np.max(sensitivity)),
            200,
        )

        trend_error = 10.0 ** (
            slope * np.log10(trend_sensitivity)
            + intercept
        )

        ax.plot(
            trend_sensitivity,
            trend_error,
            linestyle='--',
            linewidth=2,
            color='black',
        )

        intercept_sign = (
            '+'
            if intercept >= 0.0
            else '−'
        )

        equation = (
            rf'$\log_{{10}}(E) = '
            rf'{slope:.3f}\log_{{10}}(S) '
            rf'{intercept_sign} {abs(intercept):.3f}$'
        )

        annotation = (
            f'{equation}\n'
            rf'$R^2 = {r_squared:.3f}$'
        )

        ax.text(
            0.05,
            0.95,
            annotation,
            transform=ax.transAxes,
            horizontalalignment='left',
            verticalalignment='top',
            fontsize=11,
            bbox={
                'facecolor': 'white',
                'edgecolor': 'none',
                'alpha': 0.8,
            },
        )

    ax.set_xlabel(
        f'PSS sensitivity across X range '
        f'(% $\\mathit{{{metastate}}}$)'
    )
    ax.set_ylabel(
        f'Absolute PSS error '
        f'(% $\\mathit{{{metastate}}}$)'
    )

    ax.grid(
        True,
        which='major',
        linestyle='--',
        alpha=0.4,
    )

    out_path = (
        Path(out)
        / 'error_sensitivity_scatter.png'
    )

    fig.tight_layout()
    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)

    return out_path


def uvvis_pss_label(
    irradiation,
    final_pss_summary,
    sensitivity_summary,
    include_sensitivity,
    metastate,
):
    """Format a UV-Vis legend label with PSS error and optional sensitivity."""

    pss_row = (
        final_pss_summary
        .loc[
            final_pss_summary[cols.OBSERVED_IRR].eq(irradiation)
        ]
        .iloc[0]
    )

    mean_pss = float(pss_row[cols.MEAN_PSS])
    std_pss = float(pss_row[cols.PSS_STD])

    label = (
        f'{irradiation} nm: '
        f'{100.0 * mean_pss:.1f}% '
        f'± {100.0 * std_pss:.1f}%'
    )

    if include_sensitivity:
        sensitivity_values = (
            sensitivity_summary
            .loc[
                sensitivity_summary[cols.OBSERVED_IRR].eq(irradiation),
                cols.MEAN_PSS,
            ]
            .dropna()
            .to_numpy(dtype=float)
        )

        # Omit the sensitivity range if no finite sensitivity results are
        # available for this irradiation wavelength
        if sensitivity_values.size > 0:
            minimum_mean_pss = float(np.min(sensitivity_values))
            maximum_mean_pss = float(np.max(sensitivity_values))

            lower_sensitivity = max(
                mean_pss - minimum_mean_pss,
                0.0,
            )
            upper_sensitivity = max(
                maximum_mean_pss - mean_pss,
                0.0,
            )

            label += (
                ' '
                rf'$^{{+{100.0 * upper_sensitivity:.1f}\%}}'
                rf'_{{-{100.0 * lower_sensitivity:.1f}\%}}$'
            )

    label += rf' $\mathit{{{metastate}}}$'

    return label


def uvvis(
    data_df,
    irr_wls,
    final_pss_summary,
    sensitivity_summary,
    metastate_spectrum,
    uvvis_range,
    metastate,
    include_sensitivity,
    out,
):
    """Plot measured UV-Vis spectra and the extrapolated metastate spectrum.

    Each irradiated spectrum is labelled with its final mean PSS, sample
    standard deviation, and optional asymmetric quantum-yield-ratio
    sensitivity. The supplied extrapolated metastate spectrum is shown with
    the absorbance envelope obtained from mean PSS ± one sample standard
    deviation.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Validated and normalised UV-Vis data.
    irr_wls : iterable of str
        Irradiation-spectrum column labels.
    final_pss_summary : pandas.DataFrame
        Final mean PSS and standard deviation for each irradiation wavelength.
    sensitivity_summary : pandas.DataFrame
        Mean PSS results across the sampled QY-ratio range.
    metastate_spectrum : pandas.DataFrame
        Average metastate spectrum and lower and upper absorbance bounds.
    uvvis_range : tuple[float, float]
        Displayed wavelength range in nm.
    metastate : str
        Label of the metastable isomer.
    include_sensitivity : bool
        Include QY-ratio sensitivity ranges in the spectrum legend labels.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path
        Path to the saved PNG file.
    """

    wavelengths = data_df[
        cols.WAVELENGTH
    ].to_numpy(dtype=float)

    color_map = irradiation_color_map(
        irr_wls
    )

    fig, ax = plt.subplots(
        figsize=(11, 6),
    )

    ax.plot(
        wavelengths,
        data_df[cols.DARK],
        color='black',
        linewidth=2,
        label='Dark',
    )

    for irradiation in irr_wls:
        label = uvvis_pss_label(
            irradiation=irradiation,
            final_pss_summary=final_pss_summary,
            sensitivity_summary=sensitivity_summary,
            include_sensitivity=include_sensitivity,
            metastate=metastate,
        )

        ax.plot(
            wavelengths,
            data_df[irradiation],
            color=color_map[irradiation],
            linewidth=2,
            label=label,
        )

    meta_wavelengths = metastate_spectrum[
        cols.META_WAVELENGTH
    ].to_numpy(dtype=float)

    mean_meta = metastate_spectrum[
        cols.MEAN_META_ABSORBANCE
    ].to_numpy(dtype=float)

    lower_meta = metastate_spectrum[
        cols.LOWER_META_ABSORBANCE
    ].to_numpy(dtype=float)

    upper_meta = metastate_spectrum[
        cols.UPPER_META_ABSORBANCE
    ].to_numpy(dtype=float)

    ax.fill_between(
        meta_wavelengths,
        lower_meta,
        upper_meta,
        color='black',
        alpha=0.15,
        linewidth=0,
        label='PSS ± uncertainty',
    )

    ax.plot(
        meta_wavelengths,
        mean_meta,
        color='black',
        linestyle=':',
        linewidth=2.5,
        label=(
            rf'Extrapolated '
            rf'$\mathit{{{metastate}}}$'
        ),
    )

    ax.set_xlim(
        float(uvvis_range[0]),
        float(uvvis_range[1]),
    )

    ax.set_xlabel(
        'Wavelength (nm)'
    )
    ax.set_ylabel(
        'Absorbance (a.u.)'
    )

    ax.grid(
        True,
        linestyle='--',
        alpha=0.4,
    )

    # Reserve space on the right for the external legend
    fig.subplots_adjust(
        left=0.09,
        right=0.58,
        bottom=0.13,
        top=0.96,
    )

    ax.legend(
        loc='upper left',
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        fontsize=10,
    )

    out_path = Path(out) / 'uv-vis.png'

    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)

    return out_path


def pss_vs_irradiation(
    final_pss_summary,
    sensitivity_summary,
    metastate,
    include_sensitivity,
    out,
):
    """Plot final PSS against irradiation wavelength.

    Final mean PSS values are shown with their accepted-value sample standard
    deviations. If requested, wider asymmetric error bars show the range of
    filtered mean PSS values across the sampled quantum-yield ratios.

    Parameters
    ----------
    final_pss_summary : pandas.DataFrame
        Final mean PSS, sample standard deviation, and accepted-value count for
        each observed irradiation wavelength.
    sensitivity_summary : pandas.DataFrame
        Mean PSS results across the sampled quantum-yield-ratio range.
    metastate : str
        Label of the metastable isomer.
    include_sensitivity : bool
        Include asymmetric QY-ratio sensitivity error bars.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path or None
        Path to the saved PNG file, or ``None`` if no finite final PSS values
        are available.
    """

    out_path = Path(out) / 'pss_vs_irradiation.png'

    if final_pss_summary.empty:
        outlog.warning(
            'PSS-versus-irradiation plot skipped because the final PSS '
            'summary is empty.'
        )
        return None

    plot_data = final_pss_summary.copy()

    # Convert the normalised irradiation labels to numeric x coordinates
    plot_data['irradiation_nm'] = (
        plot_data[cols.OBSERVED_IRR]
        .map(utils.wavelength_value_to_float)
    )

    plot_data[cols.MEAN_PSS] = pd.to_numeric(
        plot_data[cols.MEAN_PSS],
        errors='coerce',
    )
    plot_data[cols.PSS_STD] = pd.to_numeric(
        plot_data[cols.PSS_STD],
        errors='coerce',
    )

    plot_data = (
        plot_data.loc[
            np.isfinite(plot_data[cols.MEAN_PSS])
        ]
        .sort_values('irradiation_nm')
        .reset_index(drop=True)
    )

    if plot_data.empty:
        outlog.warning(
            'PSS-versus-irradiation plot skipped because no finite final '
            'PSS values are available.'
        )
        return None

    x_values = plot_data['irradiation_nm'].to_numpy(dtype=float)
    mean_values = plot_data[cols.MEAN_PSS].to_numpy(dtype=float)
    std_values = plot_data[cols.PSS_STD].to_numpy(dtype=float)

    # A one-value result has a NaN sample standard deviation. Retain its mean
    # marker but omit its standard-deviation bar
    std_errors = np.where(
        np.isfinite(std_values) & (std_values >= 0.0),
        std_values,
        0.0,
    )

    sensitivity_errors = None

    if include_sensitivity:
        sensitivity_bounds = (
            sensitivity_summary
            .groupby(cols.OBSERVED_IRR, sort=False)[cols.MEAN_PSS]
            .agg(['min', 'max'])
            .rename(
                columns={
                    'min': '_sensitivity_min',
                    'max': '_sensitivity_max',
                }
            )
            .reset_index()
        )

        plot_data = plot_data.merge(
            sensitivity_bounds,
            on=cols.OBSERVED_IRR,
            how='left',
            sort=False,
        )

        sensitivity_min = plot_data[
            '_sensitivity_min'
        ].to_numpy(dtype=float)

        sensitivity_max = plot_data[
            '_sensitivity_max'
        ].to_numpy(dtype=float)

        lower_errors = mean_values - sensitivity_min
        upper_errors = sensitivity_max - mean_values

        # If X=1 lies slightly outside the requested sensitivity range, one
        # side can be negative. Treat that side as having zero extension from
        # the X=1 result
        lower_errors = np.where(
            np.isfinite(lower_errors) & (lower_errors > 0.0),
            lower_errors,
            0.0,
        )
        upper_errors = np.where(
            np.isfinite(upper_errors) & (upper_errors > 0.0),
            upper_errors,
            0.0,
        )

        sensitivity_errors = np.vstack(
            [lower_errors, upper_errors]
        )

    with plt.rc_context({'font.size': 15}):
        fig, ax = plt.subplots(figsize=(6, 4))

        y_axis_min = -0.1
        y_axis_max = 1.1

        ax.set_ylim(y_axis_min, y_axis_max)

        # Retain the stem-style appearance of the original plot
        ax.vlines(
            x_values,
            y_axis_min,
            mean_values,
            color='black',
            linewidth=1.2,
            linestyles=':',
            zorder=1,
        )

        ax.plot(
            x_values,
            mean_values,
            color='0.55',
            linewidth=1.2,
            zorder=2,
        )

        if (
            sensitivity_errors is not None
            and np.any(sensitivity_errors > 0.0)
        ):
            ax.errorbar(
                x_values,
                mean_values,
                yerr=sensitivity_errors,
                fmt='none',
                ecolor='tab:blue',
                capsize=8,
                elinewidth=1.2,
                alpha=0.85,
                label='QY-ratio sensitivity range',
                zorder=3,
            )

        ax.errorbar(
            x_values,
            mean_values,
            yerr=std_errors,
            fmt='.',
            color='black',
            ecolor='red',
            capsize=4,
            elinewidth=1.2,
            label='Accepted-value standard deviation',
            zorder=4,
        )

        ax.yaxis.set_major_formatter(
            PercentFormatter(xmax=1.0)
        )

        ax.set_xlabel('Irradiation wavelength (nm)')
        ax.set_ylabel(
            f'PSS (% $\\mathit{{{metastate}}}$)'
        )

        ax.set_xticks(x_values)
        ax.set_xticklabels(
            [
                utils.wavelength_label(value)
                for value in x_values
            ],
            rotation=90,
            fontsize=8 if len(x_values) > 12 else 9,
        )

        ax.grid(
            axis='y',
            linestyle='--',
            alpha=0.4,
        )

        ax.legend(fontsize=8)

        fig.tight_layout()
        fig.savefig(
            out_path,
            dpi=300,
            bbox_inches='tight',
        )
        plt.close(fig)

    return out_path


def singlepair(
    irr1,
    irr2,
    x_array,
    irr1_pss,
    irr2_pss,
    metastate,
    out,
):
    """Plot quantum-yield-ratio sensitivity for one Fischer pair.

    Two panels show the calculated PSS at irradiation wavelengths ``irr1`` and
    ``irr2`` across the sampled quantum-yield-ratio range. Non-finite results
    are omitted from each curve.

    Parameters
    ----------
    irr1, irr2 : str
        Irradiation wavelength labels defining the Fischer pair.
    x_array : array-like
        Quantum-yield-ratio values.
    irr1_pss, irr2_pss : array-like
        Calculated PSS values for irradiation wavelengths ``irr1`` and
        ``irr2`` at each quantum-yield ratio.
    metastate : str
        Label of the metastable isomer used in the y-axis labels.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path or None
        Path to the saved PNG file, or ``None`` if either panel contains no
        finite sensitivity values.
    """

    filename = (f'{irr1}_{irr2}nm_qy_ratio_sensitivity.png')
    out_path = Path(out) / filename

    plot_data = [
        {
            'target': irr1,
            'partner': irr2,
            'x': x_array,
            'y': irr1_pss,
        },
        {
            'target': irr2,
            'partner': irr1,
            'x': x_array,
            'y': irr2_pss,
        },
    ]

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(8, 8),
        sharex=False,
    )

    for ax, data in zip(axes, plot_data):
        mask = np.isfinite(data['x']) & np.isfinite(data['y'])

        if not np.any(mask):
            outlog.warning(
                f'No finite QY-sensitivity values for '
                f'{data["target"]} nm sensitivity'
            )
            plt.close(fig)
            return None
        
        ax.plot(data['x'][mask], data['y'][mask], lw=2)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))

        ax.set_xlabel(
            'Quantum-yield ratio '
            f'(X = QY ratio at {irr1} nm / QY ratio at {irr2} nm)'
        )

        ax.set_ylabel(
            f'{data["target"]} nm PSS (% $\\mathit{{{metastate}}}$)'
        )

        ax.grid(True, linestyle='--', alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    return out_path


