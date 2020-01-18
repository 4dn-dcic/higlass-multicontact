import ipywidgets as widgets
import numpy as np
import pandas as pd

from collections import OrderedDict


def get_chrom_sizes(filepath):
    df = pd.read_table(filepath, names=['chrom', 'size'], index_col='chrom')
    return df.to_dict(into=OrderedDict).get('size')


def get_selection_widgets():
    x_from = widgets.IntText(value=None, description="From:")
    x_to = widgets.IntText(value=None, description="To:")
    select_mode = widgets.ToggleButton(value=False, description="Select Mode")

    return select_mode, x_from, x_to


def link_selection_widgets(widget, select_mode, x_from, x_to):
    def handle_selection_change(selection):
        x_from.value = selection.new[0]
        x_to.value = selection.new[1]

    widget.observe(handle_selection_change, names=["selection"])

    widgets.jslink((widget, "select_mode"), (select_mode, "value"))


def lighten_color(color, amount=1):
    """
    Lightens the given color by multiplying luminosity by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 1.3) # brighten
    >> lighten_color('#F034A3', 0.6) # darken
    >> lighten_color((.3,.55,.1), 0.5) # darken
    """
    from matplotlib.colors import cnames, to_rgb
    from colorsys import rgb_to_hls, hls_to_rgb

    try:
        c = cnames[color]
    except:
        c = color
    c = rgb_to_hls(*to_rgb(c))
    return hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def get_anchor(widget_x_from, widget_x_to):
    return int(widget_x_from.value + ((widget_x_to.value - widget_x_from.value) / 2))


def rle_encode(x, dropna: bool = False):
    """Run length encoding.

    Parameters:
        x (1D array_like): sorted input array to encode
        stops (bool, optional): drop all runs of NaNs

    Returns
        numpy.ndarray: start positions

    For example:

        >>> x = [1, 1, 1, 3, 3, 5]
        >>> vrange(x)
        array([0, 3, 4])

    From:
        https://gist.github.com/nvictus/66627b580c13068589957d6ab0919e66

    """
    where = np.flatnonzero
    x = np.asarray(x)
    n = len(x)
    if n == 0:
        return np.array([], dtype=int)

    starts = np.r_[0, where(~np.isclose(x[1:], x[:-1], equal_nan=True)) + 1]

    if dropna:
        mask = ~np.isnan(x[starts])
        starts = starts[mask]

    return starts


def unique_at_freq(x, freq: int = 1):
    """Run length encoding.

    Parameters:
        x (1D array_like): sorted input array to find unique values
        freq (int): frequency at which the values should appear in x

    Returns
        numpy.ndarray: unique values at the frequency

    For example:

        >>> x = [1, 1, 1, 3, 3, 5]
        >>> unique_at_freq(x, 2)
        array([1, 3])

    """
    x = np.asarray(x)
    n = len(x)
    rle = rle_encode(x)
    values = x[rle]
    lengths = np.diff(np.r_[rle, n])

    return values[lengths >= freq]


def vrange(starts, stops, dtype = None):
    """Create concatenated ranges of integers for multiple starts/stops

    Parameters:
        starts (1D array_like): starts for each range
        stops (1D array_like): stops for each range (same shape as starts)

    Returns:
        numpy.ndarray: concatenated ranges

    For example:

        >>> starts = [1, 3, 4, 6]
        >>> stops  = [1, 5, 7, 6]
        >>> vrange(starts, stops)
        array([3, 4, 4, 5, 6])

    From:
        https://codereview.stackexchange.com/a/84980

    """
    starts = np.asarray(starts)
    stops = np.asarray(stops)
    l = stops - starts # Lengths of each range.

    if dtype:
        return np.repeat(stops - l.cumsum(), l).astype(dtype) + np.arange(l.sum()).astype(dtype)

    return np.repeat(stops - l.cumsum(), l) + np.arange(l.sum())

