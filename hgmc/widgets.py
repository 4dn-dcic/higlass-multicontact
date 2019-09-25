import higlass as hg
from IPython.display import display

from .tilesets import Hgmc1dData
from .utils import enable_selection, get_anchor

BAR_TRACK_CONFIG = {
    "uid": "bars",
    "track_type": "horizontal-bar",
    "position": "top",
    "height": 128,
    "options": {
        "colorRange": ["#ffbb33", "#e5001c", "black"],
        "labelColor": "black",
        "backgroundColor": "white",
    },
}

ANCHORS_TRACK_CONFIG = {
    "uid": "anchors",
    "track_type": "horizontal-bar",
    "position": "top",
    "height": 128,
    "options": {
        "colorRange": ["#ffbb33", "#e5001c", "black"],
        "labelColor": "black",
        "backgroundColor": "white",
    },
}


class McLevel:
    def __init__(self, tileset, higlass, x_from, x_to, select_mode, anchor=None):
        self.tileset = tileset
        self.anchor = anchor
        self.higlass = higlass
        self.x_from = x_from
        self.x_to = x_to
        self.select_mode = select_mode

    @property
    def root(self):
        return self.anchor is None


class Hgmc1d:
    def __init__(self, filepath, name: str, ignored_anchor_padding: int = 10):
        self.filepath = filepath
        self.name = name
        self.ignored_anchor_padding = ignored_anchor_padding
        self.data = Hgmc1dData(
            filepath, name=self.name, ignored_anchor_padding=self.ignored_anchor_padding
        )
        self.base_tileset = self.data.tileset()
        self.axis = hg.Track("top-axis", uid="axis")
        self.levels = []

        self.init()

    @property
    def level(self):
        return len(self.levels) - 1

    @property
    def anchors(self):
        anchors = []
        for level in self.levels:
            if level.anchor is not None:
                anchors.append(level.anchor)
        return anchors

    def get_anchor_regions(self, additional_anchors=[]):
        anchors = self.anchors + additional_anchors
        anchor_regions = []
        for anchor in anchors:
            if anchor is not None:
                anchor_regions.append(
                    [
                        anchor - self.ignored_anchor_padding,
                        anchor + self.ignored_anchor_padding + 1,
                    ]
                )
        return anchor_regions

    def show(self, level: int = None):
        if level is None:
            level = self.level

        level = max(0, min(self.level + 1, level))

        try:
            mc_level = self.levels[level]
        except IndexError:
            mc_level = self.next_level()

        display(mc_level.select_mode, mc_level.higlass, mc_level.x_from, mc_level.x_to)

    def show_next_level(self):
        self.next_level()
        self.show()

    def next_level(self):
        current_mc_level = self.levels[self.level]

        anchor = get_anchor(current_mc_level.x_from, current_mc_level.x_to)
        # start, end = anchor

        # if start is None or end is None:
        #     return current_mc_level

        self.data.add_anchors(anchor)

        tileset = self.data.tileset()

        overlays = self.get_anchor_regions([anchor])
        overlays = [
            {
                "includes": ["axis", "bars"],
                "options": {"extent": overlays, "minWidth": 6},
            }
        ]

        print(overlays)

        higlass, _, _ = hg.display(
            [
                hg.View(
                    [self.axis, hg.Track(tileset=tileset, **BAR_TRACK_CONFIG)],
                    overlays=overlays,
                )
            ]
        )

        select_mode, x_from, x_to = enable_selection(higlass)

        next_mc_level = McLevel(
            anchor=anchor,
            tileset=tileset,
            higlass=higlass,
            select_mode=select_mode,
            x_from=x_from,
            x_to=x_to,
        )

        self.levels.append(next_mc_level)

        return next_mc_level

    def init(self):
        if self.level > -1:
            return

        print

        higlass, _, _ = hg.display(
            [
                hg.View(
                    [self.axis, hg.Track(tileset=self.base_tileset, **BAR_TRACK_CONFIG)]
                )
            ]
        )

        select_mode, x_from, x_to = enable_selection(higlass)

        self.levels.append(
            McLevel(
                tileset=self.base_tileset,
                higlass=higlass,
                select_mode=select_mode,
                x_from=x_from,
                x_to=x_to,
            )
        )
