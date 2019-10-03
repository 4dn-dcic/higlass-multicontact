import ipywidgets as widgets


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
