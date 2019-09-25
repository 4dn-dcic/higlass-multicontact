import ipywidgets as widgets

def enable_selection(widget):
    x_from = widgets.IntText(value=None, description='From:')
    x_to = widgets.IntText(value=None, description='To:')

    def handle_selection_change(selection):
        x_from.value = selection.new[0]
        x_to.value = selection.new[1]

    widget.observe(handle_selection_change, names=['selection'])

    select_mode = widgets.ToggleButton(value=False, description='Select Mode')
    widgets.jslink((widget, 'select_mode'), (select_mode, 'value'))

    return select_mode, x_from, x_to

def get_anchor(widget_x_from, widget_x_to):
    return int(widget_x_from.value + ((widget_x_to.value - widget_x_from.value) / 2))
