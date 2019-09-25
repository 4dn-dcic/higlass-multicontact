import h5py
import higlass as hg
import math
import numpy as np
import warnings

# import pandas as pd

from clodius.tiles.format import format_dense_tile


def tileset(
    data,
    max_zoom,
    max_pos,
    tsinfo,
    tile_size=1024,
    ignored_regions=[],
    aggregator=np.nanmean,
    name="Some data",
    **kwargs
):
    for region in ignored_regions:
        start, end = region
        data[start:end] = np.nan

    def generate_tile(z, x):
        """
        Return tiles at the given positions.
        Parameters
        -----------
        z: int
            The zoom level (0 corresponds to most zoomed out)
        x: int
            The x tile position
        """

        tile_width = 2 ** (max_zoom - z) * tile_size

        x_start = x * tile_width
        x_end = min(max_pos, x_start + tile_width)

        # We have to convert to float to be able to use `np.nan`
        tile_data = data[x_start:x_end]

        num_to_aggregate = 2 ** (max_zoom - z)

        # add some data so that the data can be divided into squares
        divisible_x_width = num_to_aggregate * math.ceil(
            tile_data.shape[0] / num_to_aggregate
        )
        divisible_x_pad = divisible_x_width - tile_data.shape[0]

        padded_data = np.pad(
            tile_data, ((0, divisible_x_pad),), "constant", constant_values=(np.nan,)
        )

        with warnings.catch_warnings():
            # Surpress those meaningless warnings about empty slices
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out_data = aggregator(padded_data.reshape((-1, num_to_aggregate)), axis=1)

        # determine how much to pad the array
        x_pad = tile_size - out_data.shape[0]

        return np.pad(out_data, ((0, x_pad)), "constant", constant_values=(np.nan,))

    def tileset_info():
        return tsinfo

    def tiles(tile_ids):
        tiles = []

        for tile_id in tile_ids:
            # decompose the tile zoom and location
            _, z, x = tile_id.split(".")

            # generate the tile
            data = generate_tile(int(z), int(x))

            # format the tile response
            tiles.append((tile_id, format_dense_tile(data)))

        return tiles

    return hg.Tileset(tileset_info=tileset_info, tiles=tiles, name=name, **kwargs)


class Hgmc1dData:
    def __init__(
        self,
        filepath,
        name: str = "HiGlass MultiContact",
        ignored_anchor_padding: int = 3,
    ):
        self.filepath = filepath
        self.ignored_anchor_padding = ignored_anchor_padding
        self.anchors = []
        self.tile_size = 1024
        self.name = name

        with h5py.File(filepath, "r") as f:
            # print("Prepare data...")
            self.bins_to_clusters = np.zeros((len(f["clusters"]["bin"]), 2)).astype(
                np.int64
            )
            self.bins_to_clusters[:, 0] = f["clusters"]["bin"][:]
            self.bins_to_clusters[:, 1] = f["clusters"]["cluster_name"][:]

            self.min_pos = int(self.bins_to_clusters[:, 0].min())
            self.max_pos = int(self.bins_to_clusters[:, 0].max())

            self.clusters_to_bins = np.zeros((len(f["clusters"]["bin"]), 2)).astype(
                np.int64
            )
            self.clusters_to_bins[:, 0] = f["clusters"]["cluster_name"][:]
            self.clusters_to_bins[:, 1] = f["clusters"]["bin"][:]

        # Sort by bins for fast access by bins
        self.bins_to_clusters = self.bins_to_clusters[
            self.bins_to_clusters[:, 0].argsort()
        ]

        # Sort by bins for fast access by clusters
        # self.clusters_to_bins = self.clusters_to_bins[self.clusters_to_bins[:,0].argsort()]

        # Bin index
        self.bin_idx = self.get_idx(self.bins_to_clusters[:, 0])

        # Cluster index
        # self.cluster_idx = self.get_idx(self.clusters_to_bins[:,0])

        # Cluster index
        self.cluster_mask = np.zeros(
            self.bins_to_clusters[:, 1].max() + 1, dtype="bool"
        )

        # Set initial filter mask
        self.filter_masks = [self.get_filter_mask()]

        self.max_zoom = math.ceil(math.log(self.max_pos / self.tile_size) / math.log(2))
        self.max_zoom = 0 if self.max_zoom < 0 else self.max_zoom

        self.tsinfo = {
            "tile_size": self.tile_size,
            "bins_per_dimension": self.tile_size,
            "min_pos": [self.min_pos],
            "max_pos": [self.max_pos],
            "max_zoom": self.max_zoom,
            "max_width": 2 ** self.max_zoom * 1024,
        }

    def get_filter_mask(self):
        return np.ones(self.bins_to_clusters.shape[0]).astype(bool)

    def get_idx(self, arr):
        # +1 because we need to include the maximum value
        # +1 because we need to prepent a bin for the `-1`th position to
        # later on easily calculate the start and end point
        rle = np.zeros(np.max(arr) + 2)
        curr = arr[0]
        rl = 0
        for i in arr:
            if i == curr:
                rl += 1
            else:
                rle[curr + 1] = rl
                rl = 1
            curr = i
        rle[arr[-1] + 1] = rl
        return np.cumsum(rle).astype(int)

    def get_clusters_by_range(self, start, end):
        start = self.bin_idx[start]
        end = self.bin_idx[end]

        return self.bins_to_clusters[start:end, 1][self.filter_masks[-1][start:end]]

    def get_bins_by_clusters(self, start, end):
        start = self.cluster_idx[start]
        end = self.cluster_idx[end]

        return self.clusters_to_bins[start:end, 1]

    def get_anchor(self, anchor):
        try:
            start, end = anchor
        except TypeError:
            start = anchor
            end = anchor + 1
            anchor = (anchor, anchor + 1)
        except ValueError:
            start = anchor[0]
            end = anchor + 1

        return (start, end), start, end

    def add_anchors(self, *anchors, clear=False):
        for anchor in anchors:
            anchor, start, end = self.get_anchor(anchor)

            if anchor in self.anchors:
                raise ValueError("Anchor already set")

            self.anchors.append(anchor)

            clusters = self.get_clusters_by_range(start, end)

            # bins = []
            # for cluster in np.nditer(clusters):
            #     bins.append(self.get_bins_by_clusters(cluster, cluster + 1))
            # bins = np.concatenate(bins, axis=0)

            self.cluster_mask[clusters] = True
            # New Boolean mask
            self.filter_masks.append(self.cluster_mask[self.bins_to_clusters[:, 1]])

            self.cluster_mask[clusters] = False
        if clear:
            i = len(anchors)
            while i > 0:
                self.anchors.pop()
                self.filter_masks.pop()
                i -= 1

    def remove_anchors(self, *anchors):
        if len(anchors) == 0:
            try:
                self.anchors.pop()
                self.filter_masks.pop()
            except IndexError:
                pass

        for anchor in anchors:
            anchor, _, _ = self.get_anchor(anchor)
            if anchor in self.anchors:
                last_anchor = self.anchors.pop()
                self.filter_masks.pop()
                keep_anchors = []
                while last_anchor != anchor:
                    keep_anchors.append(last_anchor)
                    last_anchor = self.anchors.pop()
                    self.filter_masks.pop()

                self.add_anchors(*keep_anchors)

    def data(self):
        return self.bins_to_clusters[self.filter_masks[-1]]

    def counts_by_bin(self):
        bins = self.data()[:, 0]

        # +2 because `np.arange` and `np.histogram` exclude the last position
        return np.histogram(bins, bins=np.arange(self.min_pos, self.max_pos + 2))[0]

    def tileset(self, **kwargs):
        data = self.counts_by_bin().astype(np.float)

        ignored_regions = []
        for anchor in self.anchors:
            start, end = anchor
            start = max(0, start - self.ignored_anchor_padding)
            end = min(data.shape[0], end + self.ignored_anchor_padding)
            ignored_regions.append((start, end))

        return tileset(
            data=data,
            max_zoom=self.max_zoom,
            max_pos=self.max_pos,
            tsinfo=self.tsinfo,
            ignored_regions=ignored_regions,
            name="{}: Level {}".format(self.name, len(self.anchors)),
        )

    def anchor_tileset(self, **kwargs):
        data = np.zeros(self.max_pos - self.min_pos)

        for anchor in self.anchors:
            start, end = anchor
            start = max(0, start - self.ignored_anchor_padding)
            end = min(data.shape[0], end + self.ignored_anchor_padding)
            data[start:end] = 1

        return tileset(
            data=data,
            max_zoom=self.max_zoom,
            max_pos=self.max_pos,
            tsinfo=self.tsinfo,
            name="{}: Level {} Anchors".format(self.name, len(self.anchors)),
        )


# def mc_1d(filepath, anchors=[], **kwargs):
#     """
#     Tileset for multicontact 1D data
#     """

#     tile_size = 1024

#     def filter_data(filepath, anchors=[]):
#         with h5py.File(filepath, "r") as f:
#             data = pd.DataFrame(
#                 {"bin": f["clusters"]["bin"], "cluster": f["clusters"]["cluster_name"]}
#             )

#             min_pos = int(data["bin"].values.min())
#             max_pos = int(data["bin"].values.max())

#         counts_by_bin = np.zeros(max_pos - min_pos + 1)

#         max_cluster_name = data["cluster"].max() + 1
#         cluster_idx = np.zeros(max_cluster_name, dtype="bool")

#         for anchor in anchors:
#             clusters = data["cluster"][data["bin"] == anchor]
#             cluster_idx[clusters] = True
#             data = data.iloc[cluster_idx[data["cluster"].values]]
#             cluster_idx[clusters] = False

#         counts = data.groupby("bin").count()
#         counts_by_bin[counts.index.values] = counts["cluster"].values

#         return counts_by_bin, min_pos, max_pos

#     data, min_pos, max_pos = filter_data(filepath, anchors)

#     not_nan_data = ~np.isnan(data)

#     max_zoom = math.ceil(math.log(max_pos / tile_size) / math.log(2))
#     max_zoom = 0 if max_zoom < 0 else max_zoom

#     tsinfo = {
#         "tile_size": tile_size,
#         "bins_per_dimension": tile_size,
#         "min_pos": [min_pos],
#         "max_pos": [max_pos],
#         "max_zoom": max_zoom,
#         "max_width": 2 ** max_zoom * 1024,
#     }

#     def generate_tile(z, x):
#         """
#         Return tiles at the given positions.
#         Parameters
#         -----------
#         z: int
#             The zoom level (0 corresponds to most zoomed out)
#         x: int
#             The x tile position
#         """

#         tile_width = 2 ** (max_zoom - z) * tile_size

#         x_start = x * tile_width
#         x_end = min(max_pos, x_start + tile_width)

#         tile_data = data[x_start:x_end]
#         tile_data

#         num_to_sum = 2 ** (max_zoom - z)

#         # add some data so that the data can be divided into squares
#         divisible_x_width = num_to_sum * math.ceil(tile_data.shape[0] / num_to_sum)
#         divisible_x_pad = divisible_x_width - tile_data.shape[0]

#         padded_data = np.pad(
#             tile_data, ((0, divisible_x_pad),), "constant", constant_values=(np.nan,)
#         )

#         out_data = np.nansum(padded_data.reshape((-1, num_to_sum)), axis=1)
#         not_nan_out_data = not_nan_data[x_start:x_end]

#         # we want to calculate the means of the data points
#         na = np.pad(
#             not_nan_out_data,
#             ((0, divisible_x_pad)),
#             "constant",
#             constant_values=(np.nan,),
#         )
#         norm_out_data = np.nansum(na.reshape((-1, num_to_sum)), axis=1)
#         out_data = out_data / (norm_out_data + 1)

#         # determine how much to pad the array
#         x_pad = tile_size - out_data.shape[0]

#         return np.pad(out_data, ((0, x_pad)), "constant", constant_values=(np.nan,))

#     def tileset_info():
#         return tsinfo

#     def tiles(tile_ids):
#         tiles = []

#         for tile_id in tile_ids:
#             # decompose the tile zoom and location
#             _, z, x = tile_id.split(".")

#             # generate the tile
#             data = generate_tile(int(z), int(x))

#             # format the tile response
#             tiles.append((tile_id, format_dense_tile(data)))

#         return tiles

#     return hg.Tileset(tileset_info=tileset_info, tiles=tiles, **kwargs)
