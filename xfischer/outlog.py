"""Logging helpers for xFischer command-line runs."""

from pathlib import Path
import sys
import logging
import subprocess
from datetime import datetime
from time import perf_counter
from textwrap import TextWrapper
import numpy as np

from . import __version__
from . import utils
from . import columns as cols


logger = logging.getLogger("xfischer")

LOG_WIDTH = 120
LOG_INDENT = '    '
LOG_CONTINUATION_INDENT = '        '

_started_at = None
_output_files = {}


class _SingleBlankLineFilter(logging.Filter):
    """Suppress consecutive blank log records."""

    def __init__(self):
        super().__init__()
        self.previous_line_was_blank = True

    def filter(self, record):
        """Allow at most one blank line between non-blank records."""

        current_line_is_blank = not str(record.getMessage()).strip()
        if current_line_is_blank and self.previous_line_was_blank:
            return False

        self.previous_line_was_blank = current_line_is_blank
        return True

def setup(out):
    """Initialise file logging in the output directory."""

    global _started_at, _output_files

    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    _started_at = datetime.now()
    log_path = Path(out) / 'xfischer.log'
    _output_files = {
        str(log_path.resolve()): ('Run log', log_path),
    }

    handler = logging.FileHandler(
        log_path,
        mode='w',
        encoding='utf-8',
    )
    handler.setFormatter(logging.Formatter('%(message)s'))
    handler.addFilter(_SingleBlankLineFilter())

    logger.addHandler(handler)
    logger.propagate = False

def close():
    """Close and remove all xFischer logging handlers."""

    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)


def mainheading(title):
    """Write a main heading to the log."""

    logger.info('')
    logger.info('***************************************************************')
    logger.info(f'*{title.center(61)}*')
    logger.info('***************************************************************')
    logger.info('')


def heading(title):
    """Write a section heading to the log."""

    logger.info('')
    logger.info('      -----------------------------------------------------------      ')
    logger.info(f'     |{title.center(59)}|     ')
    logger.info('      -----------------------------------------------------------      ')
    logger.info('')

def subheading(title):
    """Write a subsection heading to the log."""

    logger.info('')
    logger.info(f'    {"-" * len(title)}')
    logger.info(f'    {title}')
    logger.info(f'    {"-" * len(title)}')
    logger.info('')


def _write_wrapped(message, level=logging.INFO, indent=LOG_INDENT):
    """Write text within the configured log width, preserving indentation."""

    for source_line in str(message).splitlines() or ['']:
        if not source_line:
            logger.log(level, '')
            continue

        leading_spaces = len(source_line) - len(source_line.lstrip(' '))
        content = source_line.lstrip(' ')
        first_indent = indent + (' ' * leading_spaces)
        continuation_indent = (
            LOG_CONTINUATION_INDENT + (' ' * leading_spaces)
        )
        wrapper = TextWrapper(
            width=LOG_WIDTH,
            initial_indent=first_indent,
            subsequent_indent=continuation_indent,
            break_long_words=False,
            break_on_hyphens=False,
        )

        for wrapped_line in wrapper.wrap(content):
            logger.log(level, wrapped_line)


def logln(message):
    """Log an indented message followed by a blank line."""
    _write_wrapped(message)
    logger.info('')


def log(message):
    """Log an indented message without a following blank line."""
    _write_wrapped(message)


def log_raw(message):
    """Log a preformatted line without automatic wrapping."""

    logger.info(LOG_INDENT + str(message))


def start_stage(title, message=None):
    """Start a timed workflow stage and return its timer value."""

    subheading(title)
    if message is not None:
        log(message)
    return perf_counter()


def finish_stage(started_at, message='Stage complete.'):
    """Log completion and elapsed wall time for a workflow stage."""

    elapsed_seconds = perf_counter() - started_at
    logln(f'{message} Elapsed time: {elapsed_seconds:.3f} seconds.')


def start_timer():
    """Return a timer value for a calculation logged on an existing line."""

    return perf_counter()


def file_saved(path, description=None):
    """Register and log a saved output file."""

    output_path = Path(path)
    output_description = description or output_path.name
    _output_files[str(output_path.resolve())] = (
        output_description,
        output_path,
    )
    logln(f'Saved {output_description}: {output_path}')


def output_summary(title='Output Summary'):
    """Log a manifest of output files written during the current run."""

    subheading(title)

    if not _output_files:
        logln('No output files were written.')
        return

    output_directory = next(
        iter(_output_files.values())
    )[1].parent.resolve()

    log(f'Output directory: {output_directory}')
    log(f'Output files written: {len(_output_files)}')

    for description, path in _output_files.values():
        if description == path.name:
            log(f'  {path.name}')
        else:
            log(f'  {description}: {path.name}')

    log('')


def intro():
    """Write the program introduction and citation information to the log."""
    logger.info('      -----------------------------------------------------------      ')
    logger.info('     |                   =====================                   |     ')
    logger.info('     |                   xFischer PSS Analysis                   |     ')
    logger.info('     |                   =====================                   |     ')
    logger.info('     |                     Hamish G. Trowell                     |     ')
    logger.info('     |                   University of Oxford                    |     ')
    logger.info('      -----------------------------------------------------------      ')
    logger.info('')
    logger.info(f'    xFischer PSS Analysis version {__version__}')
    logger.info('    Cite this work as:')
    logger.info('    * TBC')
    logger.info('')
    log(
        'Many thanks to Dominic Schatz and Davia Prischich for benchmarking '
        'this method against NMR and HPLC, and to Matthew Fuchter for '
        'supporting its development. Please see the above publication for '
        'the benchmarked results.'
    )
    logger.info('')
    logger.info('    For the original Fischer method:')
    logger.info('    E. Fischer, J. Phys. Chem., 1967, 71, 11, 3704-3706.')
    logger.info('    DOI: 10.1021/j100870a063')
    logger.info('')
    logger.info(f'    Timestamp: {datetime.now().isoformat(timespec="seconds")}')
    logger.info('')


def complete():
    """Write the successful-completion banner to the log."""

    finished_at = datetime.now()
    output_summary()

    logger.info(f'    Finished: {finished_at.isoformat(timespec="seconds")}')

    if _started_at is not None:
        duration_seconds = (finished_at - _started_at).total_seconds()
        logger.info(f'    Elapsed time: {duration_seconds:.1f} seconds')

    logger.info('')
    logger.info('***************************************************************')
    logger.info('*************** xFISCHER PSS ANALYSIS FINISHED! ***************')
    logger.info('***************************************************************')


def _explicit_cli_argument_groups(args):
    """Return grouped CLI arguments with every value-taking setting explicit."""

    argument_groups = [
        ['--input', str(args.input)],
        ['--out', str(args.out)],
    ]

    if args.singlepair:
        argument_groups.append(
            ['--singlepair', *map(str, args.singlepair)]
        )

    for flag, values in (
        ('--qyratio-range', args.qyratio_range),
        ('--darkmax-range', args.darkmax_range),
        ('--meta-range', args.meta_range),
        ('--pss-range', args.pss_range),
        ('--uvvis-range', args.uvvis_range),
    ):
        argument_groups.append([flag, *map(str, values)])

    argument_groups.extend([
        ['--meta-min', str(args.meta_min)],
        ['--hdi', str(args.hdi)],
        ['--hdi-sigma', str(args.hdi_sigma)],
    ])

    if args.violin_irrs is not None:
        violin_values = (
            [args.violin_irrs]
            if isinstance(args.violin_irrs, str)
            else args.violin_irrs
        )
        argument_groups.append(
            ['--violin-irrs', *map(str, violin_values)]
        )

    if args.sensitivity:
        argument_groups.append(['--sensitivity'])

    argument_groups.append(['--metastate', str(args.metastate)])

    return argument_groups


def _write_powershell_command(argument_groups):
    """Write a copy-pasteable PowerShell command within the log width."""

    line_indent = '      '
    continuation_indent = '        '
    available_width = LOG_WIDTH - len(continuation_indent) - 2
    lines = []
    current_line = 'python -m xfischer'

    for group in argument_groups:
        fragment = subprocess.list2cmdline(group)
        candidate = f'{current_line} {fragment}'

        if len(candidate) <= available_width:
            current_line = candidate
        else:
            lines.append(current_line)
            current_line = fragment

    lines.append(current_line)

    for index, line in enumerate(lines):
        indent = line_indent if index == 0 else continuation_indent
        continuation = ' `' if index < len(lines) - 1 else ''
        logger.info(f'{indent}{line}{continuation}')


def settings(args):
    """Log the command-line invocation and validated analysis settings.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    """

    log('Command line as entered:')
    _write_wrapped(
        utils.cli_command(),
        indent='      ',
    )
    logger.info('')

    log('Command line (verbose):')
    _write_powershell_command(_explicit_cli_argument_groups(args))
    logger.info('')

    log('Settings:')
    for key, value in vars(args).items():
        log(f'  {key}: {value}')
    logger.info('')


def errorln(message):
    """Log an error message followed by a blank line."""
    _write_wrapped(f'ERROR: {message}', level=logging.ERROR)
    logger.info('')


def error(message):
    """Log an error message without a following blank line."""
    _write_wrapped(f'ERROR: {message}', level=logging.ERROR)


def warning(message):
    """Log a warning without a following blank line."""
    _write_wrapped(f'WARNING: {message}', level=logging.WARNING)


def error_termination():
    """Write the error-termination banner to the log."""

    finished_at = datetime.now()
    output_summary('Partial Output Summary')
    logger.error(f'    Terminated: {finished_at.isoformat(timespec="seconds")}')

    if _started_at is not None:
        duration_seconds = (finished_at - _started_at).total_seconds()
        logger.error(f'    Elapsed time: {duration_seconds:.1f} seconds')

    logger.error('')
    logger.error('**********')
    logger.error('xFischer terminated with error')
    logger.error('**********')


def col_found(col_name, found_col):
    """Log detection and normalisation of a required input column.

    Parameters
    ----------
    col_name : str
        Standardised column name.
    found_col : str
        Original input column name.
    """

    if found_col and found_col != col_name:
        log(f'"{col_name}" column found from input column "{found_col}".')
    else:
        log(f'"{col_name}" column found.')


def irrs_found(irr_colnames):
    """Log the detected irradiation wavelength labels.

    Parameters
    ----------
    irr_colnames : iterable of str
        Normalised irradiation wavelength labels.
    """

    log('Irradiation wavelengths found (nm):')
    logln(f'  {str([float(i) for i in irr_colnames])}')
    

def data_overview(df, irr_wls, uvvis_range, darkmax_range):
    """Log an overview of the validated UV-Vis data.

    Parameters
    ----------
    df : pandas.DataFrame
        Validated and normalised UV-Vis data.
    irr_wls : list[str]
        Irradiation wavelength labels.
    uvvis_range : tuple[float, float]
        Wavelength range, in nm, used for absorbance-range reporting.
    darkmax_range : tuple[float, float]
        Wavelength range, in nm, used to locate the dark-spectrum absorbance
        maximum.
    """

    subheading('Data Overview')

    num_wl = len(df)
    min_wl = min(df[cols.WAVELENGTH])
    max_wl = max(df[cols.WAVELENGTH])
    spacing_wl = (max_wl - min_wl) / (num_wl - 1)
    log(f'Number of spectral wavelength points: {num_wl}')
    log(f'Wavelength range: {min_wl:.1f} - {max_wl:.1f} nm')
    logln(f'Wavelength spacing: {spacing_wl:.1f} nm')

    num_irrs = len(irr_wls)
    log(f'Number of irradiation wavelengths: {num_irrs}')
    logln(f'Observed irradiation wavelengths: {", ".join(irr_wls)} nm')

    num_irr_pairs = num_irrs * (num_irrs - 1) // 2
    num_irr_triples = num_irr_pairs * num_irrs
    log(f'Maximum available direct irradiation wavelength pairs: {num_irr_pairs}')
    logln(f'Maximum available datapoints including third-wavelength imputation: {num_irr_triples}')

    uvvis_min, uvvis_max = uvvis_range
    uvvis_mask = df[cols.WAVELENGTH].between(
        uvvis_min,
        uvvis_max,
        inclusive='both',
    )
    df_uvvis = df.loc[uvvis_mask]
    
    _, lambda_darkmax = utils.dark_lambda_max(df, darkmax_range)
    logln(
        f'lambda_max(Dark): {lambda_darkmax:.1f} nm '
        f'(search range: {darkmax_range[0]:.1f}-{darkmax_range[1]:.1f} nm)'
    )

    log(f'Absorbance ranges within UV-Vis plot range '
        f'({uvvis_min:.1f} - {uvvis_max:.1f} nm):')

    dark_min_abs = df_uvvis[cols.DARK].min()
    dark_max_abs = df_uvvis[cols.DARK].max()
    log(f'    Dark: {dark_min_abs:.2f} to {dark_max_abs:.2f} AU')

    for irr in irr_wls:
        irr_min_abs = df_uvvis[irr].min()
        irr_max_abs = df_uvvis[irr].max()
        log(f'    {irr} nm: {irr_min_abs:.2f} to {irr_max_abs:.2f} AU')
    
    log('')


def input_analysis_details(df, irr_wls, darkmax_range):
    """Log spectral values used at key wavelengths in the analysis."""

    subheading('Key Input Values')

    darkmax_index, lambda_darkmax = utils.dark_lambda_max(
        df,
        darkmax_range,
    )
    darkmax_absorbance = float(df.at[darkmax_index, cols.DARK])

    logln(
        'Dark-state reference used by the Fischer analysis: '
        f'wavelength={lambda_darkmax:.6g} nm; '
        f'A(dark)={darkmax_absorbance:.12g} AU.'
    )

    log('Spectral values at the nearest measured irradiation wavelengths:')
    widths = (14, 14, 18, 24, 18)
    headers = (
        'Requested (nm)',
        'Measured (nm)',
        'A(dark) (AU)',
        'A(irradiated) (AU)',
        'Delta A (AU)',
    )
    log_raw(
        ' | '.join(
            f'{header:>{width}}'
            for header, width in zip(headers, widths)
        )
    )
    log_raw('-+-'.join('-' * width for width in widths))

    wavelengths = df[cols.WAVELENGTH].to_numpy(dtype=float)
    for irradiation in irr_wls:
        position = utils.nearest_index(wavelengths, float(irradiation))
        row = df.iloc[position]
        measured_wavelength = float(row[cols.WAVELENGTH])
        dark_absorbance = float(row[cols.DARK])
        irradiated_absorbance = float(row[irradiation])
        delta_absorbance = irradiated_absorbance - dark_absorbance

        values = (
            str(irradiation),
            f'{measured_wavelength:.6g}',
            f'{dark_absorbance:.12g}',
            f'{irradiated_absorbance:.12g}',
            f'{delta_absorbance:.12g}',
        )
        log_raw(
            ' | '.join(
                f'{value:>{width}}'
                for value, width in zip(values, widths)
            )
        )

    log('')


def analysis_plan(
    mode,
    irradiation_wavelengths,
    pair_count,
    x_values,
):
    """Log the size and mode of the analysis about to be run."""

    x_values = np.asarray(x_values, dtype=float)
    log(f'Analysis mode: {mode}')
    log(
        f'Irradiation wavelengths ({len(irradiation_wavelengths)}): '
        f'{", ".join(str(value) for value in irradiation_wavelengths)} nm'
    )
    log(f'Unique Fischer pairs: {pair_count}')
    log(
        f'Quantum-yield ratios: {len(x_values)} values from '
        f'{x_values[0]:.12g} to {x_values[-1]:.12g}'
    )
    log('')


def _finite_text(value, scale=1.0, suffix='', precision=12):
    """Format a finite number for diagnostic logging."""

    numeric_value = float(value)
    if not np.isfinite(numeric_value):
        return 'undefined'
    return f'{numeric_value * scale:.{precision}g}{suffix}'


def _filter_text(value):
    """Format a Boolean filter result without treating NaN as true."""

    if isinstance(value, (bool, np.bool_)) and bool(value):
        return 'PASS'
    return 'FAIL'


def pss_matrix_details(pss_matrix, title, include_filters=False):
    """Log every calculated PSS row and, optionally, its filter results."""

    subheading(title)
    finite_count = int(
        np.isfinite(pss_matrix[cols.PSS].to_numpy(dtype=float)).sum()
    )
    log(
        f'Rows: {len(pss_matrix)} total; {finite_count} finite; '
        f'{len(pss_matrix) - finite_count} non-finite.'
    )

    pair_columns = [cols.PAIR_IRR1, cols.PAIR_IRR2]
    for (irr1, irr2), pair_rows in pss_matrix.groupby(
        pair_columns,
        sort=False,
    ):
        log(f'Fischer pair {irr1}/{irr2} nm:')

        for _, row in pair_rows.iterrows():
            observed = row[cols.OBSERVED_IRR]
            source = (
                'direct'
                if observed in (irr1, irr2)
                else 'third-wavelength imputation'
            )
            pss_fraction = _finite_text(row[cols.PSS])
            pss_percent = _finite_text(
                row[cols.PSS],
                scale=100.0,
                suffix='%',
            )
            details = (
                f'  observed={observed} nm ({source}); '
                f'PSS={pss_fraction} ({pss_percent})'
            )

            if include_filters:
                pss_result = _filter_text(row[cols.PSS_FILTER])
                meta_result = _filter_text(row[cols.METASTATE_FILTER])
                hdi_result = _filter_text(row[cols.HDI_FILTER])
                accepted = (
                    pss_result == 'PASS'
                    and meta_result == 'PASS'
                    and hdi_result == 'PASS'
                )
                minimum_meta = _finite_text(
                    row[cols.MIN_META_ABSORBANCE],
                    suffix=' AU',
                )
                details += (
                    f'; PSS-range={pss_result}; '
                    f'minimum metastate={minimum_meta}; '
                    f'metastate={meta_result}; HDI={hdi_result}; '
                    f'overall={"ACCEPTED" if accepted else "REJECTED"}'
                )

            if include_filters:
                # Keep one complete result per row, even when it extends past
                # the standard log width. This makes pairwise filter decisions
                # easier to scan and compare horizontally in a text editor.
                log_raw(details)
            else:
                log(details)

    log('')


def hdi_details(diagnostics):
    """Log the complete HDI calculation for each irradiation wavelength."""

    subheading('HDI Diagnostics')

    if diagnostics.empty:
        logln('No HDI diagnostics were produced.')
        return

    for _, row in diagnostics.iterrows():
        irradiation = row[cols.OBSERVED_IRR]
        log(f'Observed irradiation {irradiation} nm:')
        log(
            f'  candidates={int(row["Number of candidates"])}; '
            f'HDI fraction={_finite_text(row["HDI fraction"])}; '
            f'core size={int(row["HDI core size"])}; '
            f'sigma multiplier={_finite_text(row["HDI sigma"])}'
        )
        log(
            f'  core min={_finite_text(row["HDI core minimum"], scale=100.0, suffix="%")}; '
            f'core max={_finite_text(row["HDI core maximum"], scale=100.0, suffix="%")}; '
            f'core mean={_finite_text(row["HDI core mean"], scale=100.0, suffix="%")}; '
            f'core std={_finite_text(row["HDI core standard deviation"], scale=100.0, suffix="%")}'
        )
        lower_bound = _finite_text(
            row['HDI lower bound'],
            scale=100.0,
            suffix='%',
        )
        upper_bound = _finite_text(
            row['HDI upper bound'],
            scale=100.0,
            suffix='%',
        )
        log(
            f'  extended bounds={lower_bound} to {upper_bound}; '
            f'pass={int(row["HDI pass"])}; fail={int(row["HDI fail"])}'
        )

    log('')


def sensitivity_step(
    x_value,
    summary,
    accepted_count,
    total_count,
):
    """Log accepted results for one multipair sensitivity ratio."""

    results = []
    for _, row in summary.iterrows():
        irradiation = row[cols.OBSERVED_IRR]
        mean_text = _finite_text(
            row[cols.MEAN_PSS],
            scale=100.0,
            suffix='%',
            precision=10,
        )
        std_text = _finite_text(
            row[cols.PSS_STD],
            scale=100.0,
            suffix='%',
            precision=8,
        )
        count = int(row[cols.N_ACCEPTED])
        results.append(
            f'{irradiation} nm: mean={mean_text}, std={std_text}, n={count}'
        )

    log(
        f'X={float(x_value):.12g}; '
        f'accepted={accepted_count}/{total_count}'
    )
    for result in results:
        log(f'  {result}')


def sensitivity_overview(summary):
    """Log aggregate bounds across the complete sensitivity calculation."""

    subheading('Sensitivity Analysis Overview')
    x_values = summary[cols.QY_RATIO].drop_duplicates().to_numpy(dtype=float)
    log(
        f'Completed {len(x_values)} X values and {len(summary)} '
        'irradiation-specific summary rows.'
    )

    for irradiation, rows in summary.groupby(cols.OBSERVED_IRR, sort=False):
        finite_means = rows.loc[
            np.isfinite(rows[cols.MEAN_PSS].to_numpy(dtype=float)),
            cols.MEAN_PSS,
        ].to_numpy(dtype=float)
        accepted_counts = rows[cols.N_ACCEPTED].to_numpy(dtype=int)

        if finite_means.size:
            pss_range = (
                f'{100.0 * finite_means.min():.10g}% - '
                f'{100.0 * finite_means.max():.10g}%'
            )
        else:
            pss_range = 'undefined'

        log(
            f'{irradiation} nm: finite mean-PSS values={finite_means.size}/'
            f'{len(rows)}; mean-PSS range={pss_range}; accepted-count range='
            f'{accepted_counts.min()}-{accepted_counts.max()}'
        )

    log('')


def numeric_matrix(title, matrix, row_labels, column_labels, units):
    """Log a labelled numerical matrix in fixed-width column blocks."""

    subheading(title)
    matrix = np.asarray(matrix, dtype=float)
    log(f'Units: {units}')

    row_label_width = 13
    value_width = 14
    columns_per_block = 6

    for block_start in range(0, len(column_labels), columns_per_block):
        block_stop = block_start + columns_per_block
        block_labels = column_labels[block_start:block_stop]
        block = matrix[:, block_start:block_stop]

        log(
            'Partner-wavelength columns: '
            + ', '.join(f'{label} nm' for label in block_labels)
        )
        log_raw(
            f'{"Target (nm)":>{row_label_width}} | '
            + ' | '.join(
                f'{str(label):>{value_width}}'
                for label in block_labels
            )
        )
        log_raw(
            f'{"-" * row_label_width}-+-'
            + '-+-'.join('-' * value_width for _ in block_labels)
        )

        for label, row in zip(row_labels, block):
            values = [
                f'{value:.10g}' if np.isfinite(value) else 'undefined'
                for value in row
            ]
            log_raw(
                f'{str(label):>{row_label_width}} | '
                + ' | '.join(
                    f'{value:>{value_width}}'
                    for value in values
                )
            )

        logger.info('')

    log('')


def singlepair_curve(results):
    """Log both calculated PSS values at every single-pair X value."""

    subheading('Single-Pair Sensitivity Values')
    irr1 = results['irr1']
    irr2 = results['irr2']

    for x_value, irr1_pss, irr2_pss in zip(
        results['x_array'],
        results['irr1_pss'],
        results['irr2_pss'],
    ):
        log(
            f'X={float(x_value):.12g}; '
            f'{irr1} nm PSS={_finite_text(irr1_pss)} '
            f'({_finite_text(irr1_pss, scale=100.0, suffix="%")}); '
            f'{irr2} nm PSS={_finite_text(irr2_pss)} '
            f'({_finite_text(irr2_pss, scale=100.0, suffix="%")})'
        )

    log('')


def filter_summary(filter_counts):
    """Log pass and fail counts for the multipair filtering pipeline."""

    total_rows = filter_counts['Total rows']
    log(f'Calculated PSS rows: {total_rows}')

    for filter_name in (
        cols.PSS_FILTER,
        cols.METASTATE_FILTER,
        cols.HDI_FILTER,
    ):
        pass_count = filter_counts.get(f'{filter_name} pass', 0)
        fail_count = filter_counts.get(f'{filter_name} fail', 0)
        log(
            f'{filter_name}: {pass_count} pass, {fail_count} fail'
        )

    logln(
        'All filters: '
        f'{filter_counts["All filters pass"]} accepted, '
        f'{filter_counts["All filters fail"]} rejected'
    )


def final_pss_summary(summary):
    """Log final accepted PSS statistics for each irradiation wavelength."""

    log('Final accepted PSS values at X = 1:')

    for _, row in summary.iterrows():
        irradiation = row[cols.OBSERVED_IRR]
        mean_pss = float(row[cols.MEAN_PSS])
        std_pss = float(row[cols.PSS_STD])
        n_accepted = int(row[cols.N_ACCEPTED])

        if np.isfinite(mean_pss):
            std_text = (
                f'{100.0 * std_pss:.2f}%'
                if np.isfinite(std_pss)
                else 'undefined'
            )
            log(
                f'  {irradiation} nm: mean={100.0 * mean_pss:.2f}%, '
                f'std={std_text}, n={n_accepted}'
            )
        else:
            log(f'  {irradiation} nm: no accepted values')

    log('')


def singlepair_summary(results):
    """Log standard and sensitivity-range PSS results for one pair."""

    log('Single-pair PSS results:')

    for target in ('irr1', 'irr2'):
        irradiation = results[target]
        pss_x1 = float(results[f'{target}_pss_X1'])
        minimum = float(results[f'{target}_min'])
        maximum = float(results[f'{target}_max'])
        pss_x1_text = (
            f'{100.0 * pss_x1:.2f}%'
            if np.isfinite(pss_x1)
            else 'undefined'
        )
        log(
            f'  {irradiation} nm: PSS(X=1)={pss_x1_text}; '
            f'sensitivity range={100.0 * minimum:.2f}% to '
            f'{100.0 * maximum:.2f}%'
        )

    log('')
