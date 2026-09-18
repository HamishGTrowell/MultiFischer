"""Export tools for MultiFischer results and settings."""

from . import __version__
from . import outlog
from . import utils
from . import columns as cols

import json
from pathlib import Path
import sys
import platform
from datetime import datetime
import matplotlib
import numpy as np
import pandas as pd
import scipy
import seaborn

def settings(args):
    """Write the run settings and environment information to JSON.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.

    Returns
    -------
    pathlib.Path
        Path to the saved ``run_settings.json`` file.
    """

    command = utils.cli_command()

    settings = {
        "program": {
            "name": "MultiFischer",
            "version": __version__,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
        },
        "command": {
            "verbatim": command,
        },
        "settings": {
            **utils.json_safe(vars(args)),
            "input": str(Path(args.input).resolve()),
            "out": str(Path(args.out).resolve()),
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "executable": sys.executable,
            "working_directory": str(Path.cwd()),
            "packages": {
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "matplotlib": matplotlib.__version__,
                "scipy": scipy.__version__,
                "seaborn": seaborn.__version__,
            },
        },
    }

    path = Path(args.out) / 'run_settings.json'

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(settings, handle, indent=2)

    outlog.file_saved(path, 'run settings')
    return path


def singlepair(
    irr1,
    irr2,
    x_array,
    irr1_pss,
    irr2_pss,
    out,
):
    """Export single-pair quantum-yield-ratio sensitivity results.

    Parameters
    ----------
    irr1, irr2 : str
        Irradiation wavelength labels defining the Fischer pair.
    x_array : array-like
        Quantum-yield-ratio values.
    irr1_pss, irr2_pss : array-like
        Calculated PSS values for irradiation wavelengths ``irr1`` and
        ``irr2`` at each quantum-yield ratio.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path
        Path to the saved CSV file.
    """
    filename = f'{irr1}_{irr2}nm_qy_ratio_sensitivity.csv'
    out_path = Path(out) / filename

    df = pd.DataFrame(
        {
            'qy_ratio': x_array,
            f'{irr1}nm_pss': irr1_pss,
            f'{irr2}nm_pss': irr2_pss,
        }
    )

    df.to_csv(out_path, index=False)

    outlog.file_saved(
        out_path,
        'Single-pair QY-ratio sensitivity data',
    )

    return out_path


def meta_matrix_wide(meta_matrix):
    """Convert long-form metastate spectra into wide format.

    Parameters
    ----------
    meta_matrix : pandas.DataFrame
        Long-form metastate spectra containing Fischer-pair, wavelength, and
        metastate-absorbance columns.

    Returns
    -------
    pandas.DataFrame
        Wide-form matrix with one metastate-absorbance column per Fischer pair.
    """

    wide = meta_matrix.pivot(
        index=cols.META_WAVELENGTH,
        columns=[cols.PAIR_IRR1, cols.PAIR_IRR2],
        values=cols.META_ABSORBANCE,
    )

    wide.columns = [
        f'metastate {irr1}-{irr2} (au)'
        for irr1, irr2 in wide.columns
    ]

    wide = wide.reset_index()

    return wide


def metastate_spectra(meta_matrix, filename, out):
    """Convert and export metastate spectra in wide format.
    
    Parameters
    ----------
    meta_matrix : pandas.DataFrame
        Long-form metastate spectra containing Fischer-pair, wavelength, and
        metastate-absorbance columns. Wavelengths are in nm and absorbances are
        in arbitrary units.
    filename : str
        Output filename, including the ``.csv`` extension.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path
        Path to the saved CSV file.
    """

    df = meta_matrix_wide(meta_matrix)

    return df_to_csv(df, filename, out)


def df_to_csv(df, filename, out):
    """Export a DataFrame to CSV.

    Parameters
    ----------
    df : pandas.DataFrame
        Data to export.
    filename : str
        Output filename, including the ``.csv`` extension.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path
        Path to the saved CSV file.
    """

    out_path = Path(out) / filename
    
    df.to_csv(out_path, index=False)
    outlog.file_saved(out_path)

    return out_path


def square_matrix(
    matrix,
    irradiation_wavelengths,
    filename,
    out,
):
    """Export a directed irradiation-pair matrix to CSV.

    Parameters
    ----------
    matrix : numpy.ndarray
        Square numeric matrix.
    irradiation_wavelengths : iterable of str
        Wavelength labels defining the rows and columns.
    filename : str
        Output CSV filename.
    out : str or pathlib.Path
        Output directory.

    Returns
    -------
    pathlib.Path
        Path to the saved CSV file.
    """

    irradiation_wavelengths = list(
        irradiation_wavelengths
    )

    df = pd.DataFrame(
        matrix,
        columns=[
            f'{irr} nm'
            for irr in irradiation_wavelengths
        ],
    )

    df.insert(
        0,
        'Target irradiation wavelength (nm)',
        irradiation_wavelengths,
    )

    return df_to_csv(
        df,
        filename,
        out,
    )
