"""High-level analysis workflow."""

import numpy as np
import pandas as pd

from . import outlog
from . import io
from . import export
from . import pss
from . import plotting
from . import filters
from . import columns as cols

def run(args):
    """Run the complete xFischer workflow.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.

    Returns
    -------
    int
        Process return code. A successful run returns zero.
    """
    
    setup(args)

    data_df, irr_wls = load_data(args)

    pss_analysis(args, data_df, irr_wls)    
    outlog.complete()

    return 0


def setup(args):
    """Write introductory, settings, and environment information."""
    started_at = outlog.start_timer()
    outlog.intro()
    
    outlog.heading('Analysis Settings')
    outlog.settings(args)
    export.settings(args)
    outlog.finish_stage(
        started_at,
        'Run settings logged and exported.',
    )


def load_data(args):
    """Load input data using the validated analysis settings.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.

    Returns
    -------
    tuple[pandas.DataFrame, list[str]]
        Normalised UV-Vis data and irradiation wavelength labels.
    """
    outlog.heading('Loading Data')
    started_at = outlog.start_stage(
        'Reading Input Data',
        f'Loading input CSV file: {args.input}',
    )

    data_df, irr_wls = io.load_data(
        args.input,
        args.uvvis_range,
        args.darkmax_range,
        args.meta_range,
        args.singlepair,
        args.violin_irrs,
    )
    outlog.input_analysis_details(
        data_df,
        irr_wls,
        args.darkmax_range,
    )
    outlog.finish_stage(
        started_at,
        'Input data loaded, normalised, and validated.',
    )
    return data_df, irr_wls


def pss_analysis(args, data_df, irr_wls):
    """Dispatch to single-pair or multipair PSS analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    data_df : pandas.DataFrame
        Validated and normalised UV-Vis data.
    irr_wls : list[str]
        Irradiation wavelength labels.
    """
    outlog.mainheading('PSS Analysis')

    if args.singlepair:
        outlog.subheading('Single-Pair Analysis')
        outlog.logln(
            'Workflow branch selected: one Fischer pair across the full '
            'quantum-yield-ratio range.'
        )

        singlepair(args, data_df)

    else:
        outlog.subheading('Multi-Pair Analysis')
        outlog.logln(
            'Workflow branch selected: all unique Fischer pairs at X = 1, '
            'followed by filtered quantum-yield-ratio sensitivity analysis.'
        )

        multipair(args, data_df, irr_wls)
    

def multipair_sensitivity(
    args,
    analysis,
):
    """Calculate filtered multipair PSS statistics across the ratio range.

    At every quantum-yield ratio, new PSS and metastate matrices are
    constructed. The PSS-range, metastate, and HDI filters are independently
    reapplied before calculating the final PSS statistics.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments containing the filter settings.
    analysis : pss.FischerMultiPairAnalysis
        Initialised multipair analysis containing the pair objects and sampled
        quantum-yield ratios.

    Returns
    -------
    pandas.DataFrame
        Mean PSS, sample standard deviation, and accepted-value count for every
        quantum-yield-ratio and observed-irradiation combination.
    """

    sensitivity_summaries = []

    for X in analysis.x_array:
        pss_matrix, meta_matrix = analysis.build_matrices(X)

        filter_pipeline = (
            filters.PSSFilterPipeline(
                pss_matrix=pss_matrix,
                meta_matrix=meta_matrix,
                pss_range=args.pss_range,
                meta_range=args.meta_range,
                meta_min=args.meta_min,
                hdi=args.hdi,
                hdi_sigma=args.hdi_sigma,
            )
            .apply_pss_range_filter()
            .apply_metastate_filter()
            .apply_hdi_filter()
        )

        accepted_pss_matrix = filter_pipeline.accepted_matrix()

        x_summary = pss.summarise_accepted_pss(
            accepted_pss_matrix,
            observed_irrs=analysis.irr_wls,
        )

        x_summary.insert(
            0,
            cols.QY_RATIO,
            float(X),
        )

        outlog.sensitivity_step(
            x_value=X,
            summary=x_summary,
            accepted_count=len(accepted_pss_matrix),
            total_count=len(pss_matrix),
        )

        sensitivity_summaries.append(x_summary)

    return pd.concat(
        sensitivity_summaries,
        ignore_index=True,
    )
    

def singlepair(args, data_df):
    """Run, plot, and export single-pair sensitivity analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    data_df : pandas.DataFrame
        Validated and normalised UV–Vis data.
    """
    irr1, irr2 = args.singlepair
    analysis_started_at = outlog.start_stage(
        'Calculate Single-Pair PSS Sensitivity',
        f'Requested Fischer pair: {irr1} nm / {irr2} nm.',
    )
    outlog.log(
        f'Quantum-yield-ratio range: {args.qyratio_range[0]} to '
        f'{args.qyratio_range[1]}.'
    )
    outlog.log('Calculating both directly observed PSS curves...')

    results = pss.singlepair(
        data_df,
        irr1,
        irr2,
        darkmax_range=args.darkmax_range,
        qyratio_range=args.qyratio_range,
    )
    outlog.singlepair_summary(results)
    outlog.singlepair_curve(results)
    outlog.finish_stage(
        analysis_started_at,
        'Single-pair calculations complete.',
    )

    plot_started_at = outlog.start_stage(
        'Create Single-Pair Plot',
        'Rendering the two PSS sensitivity curves.',
    )
    plotting_singlepair_outpath = plotting.singlepair(
        irr1,
        irr2,
        results['x_array'],
        results['irr1_pss'],
        results['irr2_pss'],
        args.metastate,
        args.out,
    )
    if plotting_singlepair_outpath is not None:
        outlog.file_saved(
            plotting_singlepair_outpath,
            'Single-pair QY-ratio sensitivity plot',
        )
    outlog.finish_stage(
        plot_started_at,
        'Single-pair plotting stage complete.',
    )

    export_started_at = outlog.start_stage(
        'Export Single-Pair Data',
        f'Exporting {len(results["x_array"])} rows of sensitivity data.',
    )
    export_singlepair_outpath = export.singlepair(
        irr1,
        irr2,
        results['x_array'],
        results['irr1_pss'],
        results['irr2_pss'],
        args.out,
    )
    outlog.finish_stage(
        export_started_at,
        f'Single-pair data export complete: {export_singlepair_outpath}.',
    )

    outlog.logln('Single-pair workflow completed successfully.')

def multipair(args, data_df, irr_wls):
    """Run and export multipair PSS and sensitivity analyses.

    The standard Fischer result is calculated and filtered at ``X = 1``.
    Multipair sensitivity statistics are then independently calculated and
    filtered across the requested quantum-yield-ratio range.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    data_df : pandas.DataFrame
        Validated and normalised UV–Vis data.
    irr_wls : list[str]
        Irradiation wavelength labels.
    """
    initialise_started_at = outlog.start_stage(
        'Initialise Multi-Pair Analysis',
        'Generating the quantum-yield-ratio grid and all unique Fischer pairs.',
    )
    x_array = pss.qyratio_array(args.qyratio_range)

    analysis = pss.FischerMultiPairAnalysis(
        data_df=data_df,
        irr_wls=irr_wls,
        darkmax_range=args.darkmax_range,
        x_array=x_array,
        X=1.0,
    )
    outlog.analysis_plan(
        mode='multi-pair',
        irradiation_wavelengths=irr_wls,
        pair_count=len(analysis.pairs),
        x_values=x_array,
    )
    outlog.finish_stage(
        initialise_started_at,
        'Multi-pair analysis objects initialised.',
    )

    standard_started_at = outlog.start_stage(
        'Standard Fischer Calculation at X = 1',
        'Calculating direct and third-wavelength-imputed PSS values and '
        'extrapolated metastate spectra for every Fischer pair.',
    )
    pss_matrix, meta_matrix = analysis.build_matrices(X=1.0)
    outlog.log(
        f'Constructed PSS matrix: {len(pss_matrix)} rows x '
        f'{len(pss_matrix.columns)} columns.'
    )
    outlog.log(
        f'Constructed metastate matrix: {len(meta_matrix)} rows x '
        f'{len(meta_matrix.columns)} columns.'
    )
    outlog.pss_matrix_details(
        pss_matrix,
        'Unfiltered PSS Results at X = 1',
    )
    outlog.finish_stage(
        standard_started_at,
        'Standard X = 1 matrices complete.',
    )

    filter_started_at = outlog.start_stage(
        'Apply PSS Filters at X = 1',
        'Applying the pair-level PSS range filter, pair-level metastate '
        'filter, and irradiation-specific HDI filter in sequence.',
    )
    outlog.log(
        f'PSS acceptance range: {args.pss_range[0]} to {args.pss_range[1]}.'
    )
    outlog.log(
        f'Metastate filter: minimum absorbance {args.meta_min} AU within '
        f'{args.meta_range[0]} to {args.meta_range[1]} nm.'
    )
    outlog.log(
        f'HDI filter: core fraction {args.hdi}; '
        f'extension {args.hdi_sigma} sample standard deviations.'
    )

    filter_pipeline = (
        filters.PSSFilterPipeline(
            pss_matrix=pss_matrix,
            meta_matrix=meta_matrix,
            pss_range=args.pss_range,
            meta_range=args.meta_range,
            meta_min=args.meta_min,
            hdi=args.hdi,
            hdi_sigma=args.hdi_sigma,
        )
        .apply_pss_range_filter()
        .apply_metastate_filter()
        .apply_hdi_filter()
    )

    pss_matrix = filter_pipeline.pss_matrix
    accepted_pss_matrix = filter_pipeline.accepted_matrix()
    outlog.filter_summary(filter_pipeline.filter_counts())
    outlog.pss_matrix_details(
        pss_matrix,
        'Per-Row Filter Decisions at X = 1',
        include_filters=True,
    )
    outlog.hdi_details(filter_pipeline.hdi_diagnostics)

    # Include all observed irradiation wavelengths so that the summary retains
    # a row with n=0 when a wavelength has no accepted results.
    final_pss_summary = pss.summarise_accepted_pss(
        accepted_pss_matrix,
        observed_irrs=irr_wls,
    )
    outlog.final_pss_summary(final_pss_summary)
    outlog.finish_stage(
        filter_started_at,
        'X = 1 filtering and final PSS summary complete.',
    )

    # Later plots require one accepted finite PSS from which to extrapolate a
    # metastate spectrum. Stop here with the actual cause if every pair was
    # rejected, rather than failing later while selecting the highest PSS.
    if not np.any(
        np.isfinite(
            final_pss_summary[cols.MEAN_PSS].to_numpy(dtype=float)
        )
    ):
        raise ValueError(
            'No finite PSS values passed all filters at X = 1. Review the '
            'input spectra and filter settings.'
        )

    violin_started_at = outlog.start_stage(
        'Create PSS-Distribution Violin Plot',
        'Selecting the irradiation wavelengths to include.',
    )
    if args.violin_irrs is None:
        violin_irrs = [
            pss.best_irr(
                data_df,
                irr_wls,
                args.darkmax_range,
            )
        ]
        outlog.log(
            f'Automatic selection chose {violin_irrs[0]} nm as the '
            'irradiation wavelength with the largest dark-to-irradiated '
            'absorbance change at the dark-state maximum.'
        )
    elif args.violin_irrs == 'all':
        violin_irrs = irr_wls
        outlog.log('Configured to include every irradiation wavelength.')
    else:
        violin_irrs = args.violin_irrs
        outlog.log(
            'Configured irradiation wavelengths: '
            f'{", ".join(str(value) for value in violin_irrs)} nm.'
        )

    violin_outpath = plotting.violin(
        pss_matrix=pss_matrix,
        observed_irrs=violin_irrs,
        metastate=args.metastate,
        pss_range=args.pss_range,
        out=args.out,
    )
    if violin_outpath is not None:
        outlog.file_saved(
            violin_outpath,
            'PSS-distribution violin plot',
        )
    outlog.finish_stage(
        violin_started_at,
        'Violin-plot stage complete.',
    )

    sensitivity_started_at = outlog.start_stage(
        'Filtered Multi-Pair Sensitivity Analysis',
        f'Rebuilding and independently filtering both matrices at each of '
        f'{len(x_array)} quantum-yield ratios. Each irradiation wavelength '
        'will be logged on a separate line for every ratio.',
    )
    sensitivity_summary = multipair_sensitivity(
        args,
        analysis,
    )
    outlog.sensitivity_overview(sensitivity_summary)
    outlog.finish_stage(
        sensitivity_started_at,
        'Filtered multipair sensitivity analysis complete.',
    )

    derived_started_at = outlog.start_stage(
        'Build Derived Spectra and Pair Matrices',
        'Extrapolating the representative metastate spectrum and calculating '
        'pairwise sensitivity and absolute-error matrices.',
    )
    meta_irr, metastate_spectrum = (
        pss.highest_pss_metastate_spectrum(
            data_df=data_df,
            final_pss_summary=final_pss_summary,
        )
    )
    mean_metastate = metastate_spectrum[
        cols.MEAN_META_ABSORBANCE
    ].to_numpy(dtype=float)
    finite_metastate = mean_metastate[np.isfinite(mean_metastate)]
    outlog.log(
        f'Representative metastate spectrum is based on the highest final '
        f'mean PSS, observed at {meta_irr} nm.'
    )
    if finite_metastate.size:
        outlog.log(
            f'Representative metastate spectrum: {len(mean_metastate)} points; '
            f'{finite_metastate.size} finite; finite absorbance range '
            f'{finite_metastate.min():.12g} to '
            f'{finite_metastate.max():.12g} AU.'
        )

    matrix_irrs, sensitivity_matrix = analysis.build_sensitivity_matrix()
    error_irrs, error_matrix = analysis.build_error_matrix(
        pss_matrix=pss_matrix,
        final_pss_summary=final_pss_summary,
    )
    outlog.numeric_matrix(
        'Pair Sensitivity Matrix',
        sensitivity_matrix,
        row_labels=matrix_irrs,
        column_labels=matrix_irrs,
        units='percentage-point PSS range across the sampled X values',
    )
    outlog.numeric_matrix(
        'Pair Absolute-Error Matrix',
        error_matrix,
        row_labels=error_irrs,
        column_labels=error_irrs,
        units='absolute percentage-point difference from the final mean PSS',
    )
    outlog.finish_stage(
        derived_started_at,
        'Derived spectrum and pair matrices complete.',
    )

    plots_started_at = outlog.start_stage(
        'Create Multi-Pair Plots',
        'Rendering spectra, PSS summaries, heatmaps, and sensitivity plots.',
    )
    uvvis_outpath = plotting.uvvis(
        data_df=data_df,
        irr_wls=irr_wls,
        final_pss_summary=final_pss_summary,
        sensitivity_summary=sensitivity_summary,
        metastate_spectrum=metastate_spectrum,
        uvvis_range=args.uvvis_range,
        metastate=args.metastate,
        include_sensitivity=args.sensitivity,
        out=args.out,
    )
    if uvvis_outpath is not None:
        outlog.file_saved(uvvis_outpath, 'UV-Vis spectra plot')

    pss_irradiation_outpath = plotting.pss_vs_irradiation(
        final_pss_summary=final_pss_summary,
        sensitivity_summary=sensitivity_summary,
        metastate=args.metastate,
        include_sensitivity=args.sensitivity,
        out=args.out,
    )
    if pss_irradiation_outpath is not None:
        outlog.file_saved(
            pss_irradiation_outpath,
            'PSS-versus-irradiation plot',
        )

    sensitivity_heatmap_outpath = plotting.plot_log_heatmap(
        matrix=sensitivity_matrix,
        irradiation_wavelengths=matrix_irrs,
        colorbar_label=(
            f'PSS sensitivity across X range '
            f'(% $\\mathit{{{args.metastate}}}$)'
        ),
        filename='pair_sensitivity_heatmap.png',
        out=args.out,
    )
    if sensitivity_heatmap_outpath is not None:
        outlog.file_saved(
            sensitivity_heatmap_outpath,
            'Pair-sensitivity heatmap',
        )

    error_heatmap_outpath = plotting.plot_log_heatmap(
        matrix=error_matrix,
        irradiation_wavelengths=matrix_irrs,
        colorbar_label=(
            f'Absolute PSS error '
            f'(% $\\mathit{{{args.metastate}}}$)'
        ),
        filename='pair_error_heatmap.png',
        out=args.out,
    )
    if error_heatmap_outpath is not None:
        outlog.file_saved(
            error_heatmap_outpath,
            'Pair-error heatmap',
        )

    error_sensitivity_outpath = plotting.error_sensitivity_scatter(
        sensitivity_matrix=sensitivity_matrix,
        error_matrix=error_matrix,
        metastate=args.metastate,
        out=args.out,
    )
    if error_sensitivity_outpath is not None:
        outlog.file_saved(
            error_sensitivity_outpath,
            'Error-sensitivity scatter plot',
        )

    sensitivity_plot_outpath = plotting.multipair_sensitivity(
        sensitivity_summary=sensitivity_summary,
        metastate=args.metastate,
        out=args.out,
    )
    if sensitivity_plot_outpath is not None:
        outlog.file_saved(
            sensitivity_plot_outpath,
            'Multipair QY-ratio sensitivity plot',
        )
    outlog.finish_stage(
        plots_started_at,
        'Multi-pair plotting stage complete.',
    )

    export_started_at = outlog.start_stage(
        'Export Multi-Pair Data',
        'Writing the complete intermediate, filtered, summary, diagnostic, '
        'and pair-matrix data products.',
    )
    outlog.log(
        f'pss_matrix.csv: {len(pss_matrix)} rows x '
        f'{len(pss_matrix.columns)} columns.'
    )
    pss_matrix_outpath = export.df_to_csv(
        pss_matrix,
        'pss_matrix.csv',
        args.out,
    )

    outlog.log(
        f'metastate_spectra.csv: {len(meta_matrix)} rows x '
        f'{len(meta_matrix.columns)} columns.'
    )
    meta_matrix_outpath = export.metastate_spectra(
        meta_matrix,
        'metastate_spectra.csv',
        args.out,
    )

    outlog.log(
        f'pss_accepted_vals.csv: {len(accepted_pss_matrix)} rows x '
        f'{len(accepted_pss_matrix.columns)} columns.'
    )
    accepted_pss_matrix_outpath = export.df_to_csv(
        accepted_pss_matrix,
        'pss_accepted_vals.csv',
        args.out,
    )

    outlog.log(
        f'hdi_results.csv: {len(filter_pipeline.hdi_diagnostics)} rows x '
        f'{len(filter_pipeline.hdi_diagnostics.columns)} columns.'
    )
    hdi_results_outpath = export.df_to_csv(
        filter_pipeline.hdi_diagnostics,
        'hdi_results.csv',
        args.out,
    )

    outlog.log(
        f'pss_final_summary.csv: {len(final_pss_summary)} rows x '
        f'{len(final_pss_summary.columns)} columns.'
    )
    pss_final_summary_outpath = export.df_to_csv(
        final_pss_summary,
        'pss_final_summary.csv',
        args.out,
    )

    outlog.log(
        f'pss_sensitivity_summary.csv: {len(sensitivity_summary)} rows x '
        f'{len(sensitivity_summary.columns)} columns.'
    )
    pss_sensitivity_summary_outpath = export.df_to_csv(
        sensitivity_summary,
        'pss_sensitivity_summary.csv',
        args.out,
    )

    outlog.log(
        f'pair_sensitivity_matrix.csv: {sensitivity_matrix.shape[0]} rows x '
        f'{sensitivity_matrix.shape[1]} columns.'
    )
    sensitivity_matrix_outpath = export.square_matrix(
        sensitivity_matrix,
        matrix_irrs,
        'pair_sensitivity_matrix.csv',
        args.out,
    )

    outlog.log(
        f'pair_error_matrix.csv: {error_matrix.shape[0]} rows x '
        f'{error_matrix.shape[1]} columns.'
    )
    error_matrix_outpath = export.square_matrix(
        error_matrix,
        matrix_irrs,
        'pair_error_matrix.csv',
        args.out,
    )
    outlog.finish_stage(
        export_started_at,
        'Multi-pair data export complete.',
    )
