"""General command-line, wavelength, and numerical utility functions."""

import subprocess
import sys
import numpy as np
from pathlib import Path

from . import columns as cols

def cli_command():
    """Return a copy-pasteable command representing the current invocation."""

    program_name = Path(sys.argv[0]).name

    if program_name == '__main__.py':
        command = [
            'python',
            '-m',
            'multifischer',
            *sys.argv[1:],
        ]
    else:
        command = [
            program_name,
            *sys.argv[1:],
        ]

    return subprocess.list2cmdline(command)

def wavelength_value_to_float(value: object):
    """Convert a wavelength-like value to a float.

    Parameters
    ----------
    value : object
        Numeric value or text containing a wavelength with an optional ``nm``
        suffix.

    Returns
    -------
    float
        Numeric wavelength.

    Raises
    ------
    ValueError
        If the value cannot be interpreted as a wavelength.
    """

    text = str(value).lower().replace('nm', '').strip()
    return float(text)


def wavelength_label(value):
    """Normalise a wavelength-like value to a compact string label.

    Integer-valued wavelengths are returned without a decimal part, so ``365``
    and ``365.0`` both become ``"365"``.

    Parameters
    ----------
    value : object
        Wavelength-like value.

    Returns
    -------
    str
        Normalised wavelength label.
    """

    wl = wavelength_value_to_float(value)
    if wl.is_integer():
        return str(int(wl))
    return f'{wl:g}'


def json_safe(value):
    """Recursively convert values into JSON-serialisable Python objects."""
    if isinstance(value, dict):
        return {key: json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def nearest_index(values, target):
    """Return the position of the finite value nearest to a target.

    Parameters
    ----------
    values : array-like
        Numeric values to search.
    target : float
        Target value.

    Returns
    -------
    int
        Zero-based position of the nearest value.
    """

    array = np.asarray(values, dtype=float)
    return int(np.nanargmin(np.abs(array - float(target))))


def dark_lambda_max(df, darkmax_range):
    """Locate the dark-spectrum absorbance maximum within a wavelength range.

    Parameters
    ----------
    df : pandas.DataFrame
        UV–Vis data containing wavelength and dark-spectrum columns.
    darkmax_range : tuple[float, float]
        Inclusive wavelength range, in nm, used for the search.

    Returns
    -------
    tuple[int, float]
        DataFrame index and corresponding wavelength of the maximum dark-state
        absorbance.
    """

    wl_data = df[cols.WAVELENGTH]

    darkmax_band = (
        (wl_data >= darkmax_range[0]) &
        (wl_data <= darkmax_range[1])
    )

    idx_darkmax = df.loc[darkmax_band, cols.DARK].idxmax()
    lambda_darkmax = float(df.at[idx_darkmax, cols.WAVELENGTH])

    return idx_darkmax, lambda_darkmax
