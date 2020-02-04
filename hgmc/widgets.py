import higlass as hg
from IPython.display import display
from matplotlib.cm import get_cmap

from .tilesets import Hgmc1dData
from .utils import get_selection_widgets, link_selection_widgets, get_anchor

cmap = get_cmap("cool")

BAR_TRACK_CONFIG = {
    "uid": "bars",
    "track_type": "horizontal-bar",
    "position": "top",
    "height": 128,
    "options": {
        "colorRange": ["#ffffe0", "#0000bf"],
        "labelColor": "black",
        "backgroundColor": "white",
    },
}

ANCHORS_OPTIONS = {
    "fillColor": "orange",
    "fillOpacity": 1,
    "outline": "white",
    "outlineWidth": 1,
    "outlinePos": ["left", "right"],
}

NEW_ANCHORS_OPTIONS = {
    "fillOpacity": 0,
    "stroke": "orange",
    "strokeWidth": 2,
    "strokePos": ["top", "left", "right"],
    "outline": "white",
    "outlineWidth": 1,
    "outlinePos": ["left", "right"],
}


def get_higlass_widget(mc_level):
    higlass, _, _ = hg.display(mc_level.higlass_views, no_fuse=True)

    link_selection_widgets(
        higlass, mc_level.select_mode, mc_level.x_from, mc_level.x_to
    )

    return higlass


class McLevel:
    def __init__(self, tileset, higlass_views, x_from, x_to, select_mode, anchor=None):
        self.tileset = tileset
        self.anchor = anchor
        self.higlass_views = higlass_views
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

        higlass_views = [
            hg.View(
                [self.axis, hg.Track(tileset=self.base_tileset, **BAR_TRACK_CONFIG)]
            )
        ]

        select_mode, x_from, x_to = get_selection_widgets()

        self.levels.append(
            McLevel(
                tileset=self.base_tileset,
                higlass_views=higlass_views,
                select_mode=select_mode,
                x_from=x_from,
                x_to=x_to,
            )
        )

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

    def show_all_levels(self, track_height: int = 36):
        tracks = [self.axis]
        overlays = []
        curr_anchors = []

        for index, level in enumerate(self.levels):

            uid = "bars-{}".format(index)

            tracks.append(
                hg.Track(
                    tileset=level.tileset,
                    **{**BAR_TRACK_CONFIG, **{"uid": uid, "height": track_height}}
                )
            )

            if level.anchor is not None:
                new_anchor = [
                    level.anchor - self.ignored_anchor_padding,
                    level.anchor + self.ignored_anchor_padding + 1,
                ]
                overlays.append(
                    {
                        "uid": "overlays-{}-new".format(index),
                        "includes": ["bars-{}".format(index - 1)],
                        "options": {**{"extent": [new_anchor]}, **NEW_ANCHORS_OPTIONS},
                    }
                )
                curr_anchors.append(new_anchor)
                if len(curr_anchors):
                    overlays.append(
                        {
                            "uid": "overlays-{}-prev".format(index),
                            "includes": [uid],
                            "options": {
                                **{"extent": curr_anchors.copy()},
                                **ANCHORS_OPTIONS,
                            },
                        }
                    )
            else:
                overlays.append(
                    {
                        "uid": "overlays-{}".format(index),
                        "includes": ["axis"],
                        "options": {
                            **{
                                "extent": self.get_anchor_regions([level.anchor]),
                                "minWidth": 4,
                            },
                            **ANCHORS_OPTIONS,
                        },
                    }
                )

        higlass, _, _ = hg.display([hg.View(tracks, overlays=overlays)], no_fuse=True)

        display(higlass)

    def show_current_level(self):
        level = self.level

        mc_level = self.levels[level]

        higlass = get_higlass_widget(mc_level)

        display(mc_level.select_mode, higlass, mc_level.x_from, mc_level.x_to)

    def show(self, level: int = None, all: bool = False):
        if all:
            self.show_all_levels()
            return

        if level is None:
            self.show_current_level()
            return

        level = max(0, min(self.level + 1, level))

        if level > self.level:
            self.show_next_level()

        mc_level = self.levels[level]

        higlass = get_higlass_widget(mc_level)

        display(mc_level.select_mode, higlass, mc_level.x_from, mc_level.x_to)

    def show_next_level(self):
        self.next_level()
        self.show_current_level()

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
                "options": {**{"extent": overlays, "minWidth": 4}, **ANCHORS_OPTIONS},
            }
        ]

        higlass_views = [
            hg.View(
                [self.axis, hg.Track(tileset=tileset, **BAR_TRACK_CONFIG)],
                overlays=overlays,
            )
        ]

        select_mode, x_from, x_to = get_selection_widgets()

        next_mc_level = McLevel(
            anchor=anchor,
            tileset=tileset,
            higlass_views=higlass_views,
            select_mode=select_mode,
            x_from=x_from,
            x_to=x_to,
        )

        self.levels.append(next_mc_level)

        return next_mc_level

    def reset(self):
        while len(self.levels) > 1:
            self.levels.pop()

        self.data.remove_anchors()
