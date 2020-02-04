import gzip
import h5py
import math
import numpy as np
import os
import time

from .utils import get_chrom_sizes, rle_encode, vrange


def tsv_to_h5(
    tsv_filepath: str,
    chrom_sizes_filepath: str,
    chrom: str,
    bin_size: int,
    cluster_size: int = 3,
    header: bool = False,
    timeit: bool = False,
    verbose: bool = False
):
    t0 = time.time()

    base, ext = os.path.splitext(tsv_filepath)

    open_tsv = lambda f: open(f, 'r')

    if ext == '.gz' or ext == '.gzip':
        open_tsv = lambda f: gzip.open(f, 'rt')
        base, ext = os.path.splitext(base)

    h5_filepath = base + ".h5"

    chrom_size = get_chrom_sizes(chrom_sizes_filepath).get(chrom)
    num_bins = math.ceil(chrom_size / bin_size)

    with h5py.File(h5_filepath, 'w') as h5:
        with open_tsv(tsv_filepath) as tsv:
            if header:
                tsv.readline() # Skip header

            # Count lines
            num_lines = 0
            for line in tsv:
                num_lines += 1

            if verbose:
                print(f'File contains {num_lines} clusters')

            cluster_to_bin = np.zeros((num_lines * cluster_size, 2)).astype(np.uint32)
            cluster_to_bin[:, 0] = np.repeat(np.arange(num_lines), cluster_size)
            cluster_offset = np.arange(0, (num_lines + 1) * cluster_size, cluster_size).astype(np.uint32)
            cluster_props = np.zeros((num_lines, 2)).astype(np.uint32) # frequency, span

            tsv.seek(0)

            if header:
                tsv.readline() # Skip header

            i = 0
            for line in tsv:
                # Expects:
                # chr1, pos1, chr2, pos2, chr3, pos3, ..., freq
                cols = line.strip().split('\t')

                # Add bins
                for j in range(cluster_size):
                    cluster_to_bin[i * cluster_size + j, 1] = np.uint32(cols[j * 2 + 1]) / bin_size

                # Add props
                cluster_props[i, 0] = np.uint32(cols[cluster_size * 2]) # frequency
                cluster_props[i, 1] = np.uint32(cols[cluster_size * 2 - 1]) - np.uint32(cols[1]) # span

                i += 1

            h5.create_dataset("cluster_to_bin", data=cluster_to_bin)
            h5.create_dataset("cluster_offset", data=cluster_offset)
            h5.create_dataset("cluster_props", data=cluster_props)

            # Sort by bins
            bin_to_cluster = np.copy(cluster_to_bin)
            bin_to_cluster[:, [0, 1]] = bin_to_cluster[:, [1, 0]]
            bin_to_cluster = bin_to_cluster[np.argsort(bin_to_cluster[:,0])]

            # Determine bin offsets.
            bin_offset = np.zeros(num_bins + 1).astype(np.uint32)

            # First we determine the run-length encoding of the bin list
            bin_rle = rle_encode(bin_to_cluster[:, 0])

            rle_bins = bin_to_cluster[:, 0][bin_rle]
            rle_bins = np.r_[bin_to_cluster[:, 0][bin_rle], num_bins]
            bin_rle = np.r_[bin_rle, bin_to_cluster[:, 0].size]

            bin_offset[rle_bins] = bin_rle

            # Some bins might have not been in any cluster. We need to handle them manually
            bin_mask = np.ones(num_bins + 1).astype(bool)
            bin_mask[rle_bins] = False
            bin_mask[num_bins] = False

            i = 0
            for untouched_bin in np.where(bin_mask)[0]:
                # For bins that are not part of any cluster, i.e., untouched we will take the offset
                # from the next bin higher that is in a cluster
                while untouched_bin >= rle_bins[i]:
                    i += 1
                bin_offset[untouched_bin] = bin_offset[rle_bins[i]]

            h5.create_dataset("bin_to_cluster", data=bin_to_cluster)
            h5.create_dataset("bin_offset", data=bin_offset)

    if timeit:
        print(f'Took {((time.time() - t0) / 60):.2f} minutes')


def verify_queried_clusters(
    h5, query: list, clusters, timeit: bool = False, verbose: bool = False
):
    t0 = time.time()

    num_total_clusters = h5['cluster_props'].shape[0]

    cluster_mask = np.zeros(num_total_clusters + 1).astype(bool)

    t1 = time.time()
    if timeit:
        print(f'Mask generation took {t1-t0:.1f} sec')

    cluster_mask[clusters] = True

    cluster_starts = h5['cluster_offset'][cluster_mask]
    cluster_stops = h5['cluster_offset'][np.roll(cluster_mask, 1)]
    cluster_lengths = cluster_stops - cluster_starts

    t2 = time.time()
    if timeit:
        print(f'Cluster start/stop extraction took {t2-t1:.1f} sec')

    ranges = vrange(cluster_starts, cluster_stops, np.uint64)

    t3 = time.time()
    if timeit:
        print(f'Vrange calculation took {t3-t2:.1f} sec')

    cluster_to_bin = h5['cluster_to_bin'][:, 1]

    t4 = time.time()
    if timeit:
        print(f'Cluster-to-bin extraction took {t4-t3:.1f} sec')

    bins = cluster_to_bin[ranges]

    t5 = time.time()
    if timeit:
        print(f'Bin extraction took {t5-t4:.1f} sec')
        print(f'Total took {t5-t0:.1f} sec')

    for i, q in enumerate(query):
        feature_mask, feature_freq = q

        cluster_features = feature_mask[bins]

        j = 0
        for length in cluster_lengths:
            cluster_feat_freq = cluster_features[j:j+length].sum()
            if cluster_feat_freq < feature_freq:
                if verbose:
                    print(
                        f'Damn! Cluster #{j} has only {cluster_feat_freq} '
                        f'bins with feature {i} but >= {feature_freq} are '
                        f'required'
                    )
                return False

            j += length

    if verbose:
        print('Hooray! The clusters conform to the query.')

    return True


def query_by_features(
    h5,
    query: list = [],
    timeit: bool = False,
    verify: bool = False,
    verbose: bool = False
):
    t0 = time.time()

    ## INIT
    num_bins = h5['bin_offset'].size - 1
    num_clusters = h5['cluster_props'].shape[0]

    ## BIN MASKS:
    # This mask will be used to determine the currently selected bins
    bin_mask = np.ones(num_bins + 1).astype(bool)
    bin_mask[-1] = False # Must always be false!
    # This mask will be used to save which bins we have looked at already
    bin_mask_visited = np.zeros(num_bins + 1).astype(bool)

    ## CLUSTER MASK
    cluster_mask = np.ones(num_clusters + 1).astype(bool)
    cluster_mask[-1] = False # Must always be false!
    cluster_mask_current = np.zeros(num_clusters + 1).astype(np.uint32)
    t1 = time.time()
    if timeit:
        print(f'Mask creation took {(t1 - t0):.2f} sec')

    bin_to_cluster = h5['bin_to_cluster'][:, 1]
    cluster_to_bin = h5['cluster_to_bin'][:, 1]
    t2 = time.time()
    if timeit:
        print(f'bin_to_cluster and cluster_to_bin extraction took {(t2 - t1):.2f} sec')

    ## RUN QUERY
    while len(query):
        t2 = time.time()
        # 1. Find bins that match the query
        feature_mask, feature_freq = query.pop()
        bin_mask[:-1] *= feature_mask

        # 2. Update `bin_mask_visited`
        # Note that `+=` is the same as `np.logical_or()`
        bin_mask_visited += bin_mask

        # 3. Get cluster offsets with contact points at `bin_mask`
        bin_starts = h5['bin_offset'][bin_mask]
        bin_stops = h5['bin_offset'][np.roll(bin_mask, 1)]
        t3 = time.time()
        if timeit:
            print(f'Bin starts and ends extraction took {(t3-t2):.2f} sec')

        # 3.1 Get cluster IDs
        cluster_mask_current[cluster_mask_current > 0] = 0 # Reset mask
        for i, start in enumerate(bin_starts):
            stop = bin_stops[i]
            clusters_at_bin = bin_to_cluster[start:stop]

            # Since this is the first round we do `+= 1`.
            cluster_mask_current[clusters_at_bin] += 1

        # Update the overall cluster mask
        # Note that `*=` is equivalent to `np.logical_and`
        cluster_mask *= cluster_mask_current >= feature_freq

        t4 = time.time()
        if timeit:
            print(f'Getting the unique cluster IDs took {(t4-t3):.2f} sec')

        # Check if this was the last query condition
        if not len(query):
            # We can already leave the loop
            break

        # 4. Get all bins from the matching clusters
        bin_mask[bin_mask] = False # Unset previous bin mask
        cluster_starts = h5['cluster_offset'][cluster_mask]
        cluster_stops = h5['cluster_offset'][np.roll(cluster_mask, 1)]

        t5 = time.time()
        if timeit:
            print(f'Cluster starts and ends extraction took {(t5-t4):.2f} sec')

        ranges = vrange(cluster_starts, cluster_stops, np.uint64)

        t6 = time.time()
        if timeit:
            print(f'Vrange calculation took {(t6-t5):.2f} sec')

        bins_in_cluster = cluster_to_bin[ranges]
        bin_mask[bins_in_cluster] = True

        t7 = time.time()
        if timeit:
            print(f'Bin extraction took {(t7-t6):.2f} sec')

        # 5. Exclude previously visited bins using XOR
        bin_mask = np.logical_xor(bin_mask, bin_mask_visited)

        t8 = time.time()
        if timeit:
            print(f'Excluding previously visited bins took {(t8-t7):.2f} sec')

    if timeit:
        print(f'Total query took {(time.time() - t0):.2f} sec')

    cluster_ids = np.where(cluster_mask)[0]

    if verify:
        if verbose:
            print('Verify results...')

        passed = verify_queried_clusters(
            h5, query, cluster_ids, timeit=timeit, verbose=verbose
        )

        assert passed, 'Results are invalid'

    return cluster_ids


def clusters_to_bins(h5, clusters, timeit: bool = False):
    t0 = time.time()

    num_clusters = h5['cluster_props'].shape[0]

    cluster_mask = clusters

    if clusters.dtype != np.dtype('bool'):
        cluster_mask = np.zeros(num_clusters + 1).astype(bool)

        t1 = time.time()
        if timeit:
            print(f'Mask generation took {t1-t0:.1f} sec')

        cluster_mask[clusters] = True

    cluster_starts = h5['cluster_offset'][cluster_mask]
    cluster_stops = h5['cluster_offset'][np.roll(cluster_mask, 1)]

    t2 = time.time()
    if timeit:
        print(f'Cluster start/stop extraction took {t2-t1:.1f} sec')

    ranges = vrange(cluster_starts, cluster_stops, np.uint64)

    t3 = time.time()
    if timeit:
        print(f'Vrange calculation took {t3-t2:.1f} sec')

    cluster_to_bin = h5['cluster_to_bin'][:]

    t4 = time.time()
    if timeit:
        print(f'Cluster-to-bin extraction took {t4-t3:.1f} sec')

    bins = cluster_to_bin[ranges]

    t5 = time.time()
    if timeit:
        print(f'Bin extraction took {t5-t4:.1f} sec')
        print(f'Total took {t5-t0:.1f} sec')

    return bins


def bins_to_cluster_mask(h5, bins):
    num_bins = h5['bin_offset'].size - 1
    num_clusters = h5['cluster_props'].shape[0]

    bin_mask = np.zeros(num_bins + 1).astype(bool)
    bin_mask[bins] = True
    cluster_mask = np.zeros(num_clusters + 1).astype(bool)

    bin_starts = h5['bin_offset'][bin_mask]
    bin_stops = h5['bin_offset'][np.roll(bin_mask, 1)]

    bin_to_cluster = h5['bin_to_cluster'][:, 1]

    for i, start in enumerate(bin_starts):
        stop = bin_stops[i]
        cluster_mask[bin_to_cluster[start:stop]] = True

    return cluster_mask

def bins_to_clusters(*args):
    return np.where(bins_to_cluster_mask(*args))[0]


def clusters_to_unique_bins(*args, **kwargs):
    return np.unique(clusters_to_bins(*args, **kwargs)[:,1])


# def feature_support(h5, features):
#     # 1. Feature > bin mask
#     # 2. bin mask > clusters
#     # 3. clusters > frequency

def support(h5, query):
    clusters = query_by_features(h5, query)
    cluster_freqs = h5['cluster_props'][:, 0]

    N = cluster_freqs.sum()
    count = cluster_freqs.size * cluster_freqs[clusters].sum()

    return (count / N, count)


def cross_support(support_count):
    return support_count.min() / support_count.max()


def cluster_cross_support(h5, clusters):
    cluster_freqs = h5['cluster_props'][:, 0][clusters]
    return cross_support(cluster_freqs)


def query_cross_support(h5, query, feature_freq):
    support_count = np.zeros(len(query)).astype(np.uint32)

    for i, q in enumerate(query):
        count, _ = support(h5, q)
        support_count[i] = count

    return cross_support(support_count)


def lift(h5):
    return None


def confidence():
    return None


def shannon_entropy(rel_feature_counts):
    log2_rel_counts = np.log2(
        rel_feature_counts,
        out=np.zeros_like(rel_feature_counts),
        where=rel_feature_counts!=0
    )
    return -(rel_feature_counts * log2_rel_counts).sum(axis=0)


def information_content(rel_feature_counts, cluster_counts):
    # Number of feature states
    s, _ = rel_feature_counts.shape

    # Small sample correction
    e = 1 / np.log(2) * np.divide(
        s - 1,
        2 * cluster_counts,
        out=np.zeros_like(cluster_counts),
        where=cluster_counts!=0
    ).mean()

    return np.log2(s) - (shannon_entropy(rel_feature_counts) + e)


def contact_feature_logo(h5, feature, verbose: bool = False):
    t0 = time.time()

    num_states, num_bins = feature.shape

    bin_offset = h5['bin_offset'][:]
    cluster_offset = h5['cluster_offset'][:]
    bin_to_cluster = h5['bin_to_cluster'][:, 1]
    cluster_to_bin = h5['cluster_to_bin'][:]
    cluster_freq = h5['cluster_props'][:, 0]

    num_bins = bin_offset.size - 1
    num_clusters = cluster_offset.size - 1

    bin_counts = np.zeros(num_bins).astype(np.uint32)
    cluster_mask = np.zeros(num_clusters + 1).astype(bool)

    feature_counts = np.zeros((num_states, num_bins))
    num_clusters_per_bin = np.zeros(num_bins)

    if verbose:
        print(f'Initialization took {(time.time() - t0):.2f} sec')

    steps = 100
    dt = np.zeros(steps)
    for i, bin in enumerate(np.arange(num_bins)):
        t1 = time.time()

        # 0. Reset cluster mask
        cluster_mask[cluster_mask] = False

        # Get the states for bins that are in contact with this bin
        # 1. Get the clusters at this bin
        bin_start = bin_offset[bin]
        bin_stop = bin_offset[bin + 1]
        num_clusters_per_bin[bin] = bin_stop - bin_start

        cluster_mask[bin_to_cluster[bin_start:bin_stop]] = True

        # 2. Get cluster bins
        cluster_starts = cluster_offset[cluster_mask]
        cluster_stops = cluster_offset[np.roll(cluster_mask, 1)]

        ranges = vrange(cluster_starts, cluster_stops, np.uint64)

        cluster_bins = cluster_to_bin[ranges]

        # Unset previous counts
        bin_counts[bin_counts > 0] = 0

        # 3. Count appearance of cluster bins
        for cb in cluster_bins:
            cluster_id, cluster_bin = cb
            bin_counts[cluster_bin] += cluster_freq[cluster_id]

        # # 3. Count appearance of cluster bins
        # unique_bins, bin_counts = np.unique(cluster_bins, return_counts=True)

        # Remove the current bin as we only want to count other bins that are
        # in contact with the current bin
        # bin_counts[unique_bins == bin] = 0
        bin_counts[bin] = 0

        # 4. Get the feature counts across the cluster bins
        feature_counts[:, bin] = (feature * bin_counts).sum(axis=1)

        dt[i % steps] = time.time() - t1

        if verbose and i % steps == steps - 1:
            print(f'The last {steps} steps took {dt.sum():.1f} sec')

    total_feature_counts = feature_counts.sum(axis=0)
    rel_feature_counts = np.divide(
        feature_counts,
        total_feature_counts,
        out=np.zeros_like(feature_counts),
        where=total_feature_counts != 0
    )

    r = information_content(rel_feature_counts, num_clusters_per_bin)

    if verbose:
        print(f'Total logo computation took {((time.time() - t0) / 60):.2f} min')

    return rel_feature_counts * r
