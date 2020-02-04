import numpy as np

from matplotlib import pyplot as plt

from .clusters import bins_to_cluster_mask, clusters_to_bins

def plot_distribution(
    Y,
    labels = None,
    Y_comp = None,
    colors = None,
    figsize: tuple = (12, 8),
    title: str = None,
    ylim: tuple = None,
    yticks = None,
    primary_label: str = None,
    comparison_label: str = None
):
    Y = np.asarray(Y)
    num_features = Y.size

    X = np.arange(num_features)

    width_primary = 4/6
    width_secondary = 1/6

    color_primary = '#00abff'
    color_secondary = '#000000'

    has_comp = Y_comp is not None

    fig, ax = plt.subplots(
        nrows=2 if has_comp else 1,
        ncols=1,
        figsize=figsize,
        sharex=True,
        gridspec_kw=dict(
            height_ratios=([1] + [0.5] * has_comp),
            hspace=0.25
        )
    )

    curr_ax = ax[0] if Y_comp is not None else ax

    curr_ax.bar(
        X,
        Y,
        width=width_primary,
        color=colors if colors is not None else color_primary,
        label=primary_label if primary_label is not None else 'Found cluster bins'
    )
    if Y_comp is not None:
        Y_comp = np.asarray(Y_comp)
        curr_ax.bar(
            X + (width_primary / 2) - (width_secondary / 2),
            Y_comp,
            width=width_secondary,
            color=color_secondary,
            alpha=0.33,
            label=comparison_label if comparison_label is not None else 'Overlapping cluster bins'
        )
        curr_ax.legend()

    if title is not None:
        curr_ax.set_title(title)
    else:
        curr_ax.set_title('Absolute Feature Counts by Clusters by Bins')
    curr_ax.set_xticks(X)
    if labels is not None:
        curr_ax.set_xticklabels(labels)
    if ylim is not None:
        curr_ax.set_ylim(*ylim)
    if yticks is not None:
        curr_ax.set_yticks(yticks)

    curr_ax.spines['top'].set_visible(False)
    curr_ax.spines['left'].set_visible(False)
    curr_ax.spines['right'].set_visible(False)

    if Y_comp is not None:
        Y_comp = np.asarray(Y_comp)
        ax[1].bar(
            X,
            np.divide(Y, Y_comp, out=np.zeros_like(Y).astype(float), where=Y_comp!=0),
            width=width_primary,
            color=colors if colors is not None else color_primary,
        )

        ax[1].set_title('Relative Feature Coverage')
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['left'].set_visible(False)
        ax[1].spines['right'].set_visible(False)


def plot_cluster_feature_distribution(h5, clusters_bins, features, figsize: tuple = (12, 8)):
    num_features = len(features)
    counts = np.zeros(num_features).astype(np.uint32)
    counts_all_clusters = np.zeros(num_features).astype(np.uint32)
    labels = []

    bins = clusters_bins[:, 1]
    clusters = clusters_bins[:, 0]

    all_clusters_bins = clusters_to_bins(
        h5,
        bins_to_cluster_mask(h5, np.unique(bins))
    )
    all_bins = all_clusters_bins[:, 1]
    all_clusters = all_clusters_bins[:, 0]

    cluster_freqs = h5['cluster_props'][:, 0][clusters]
    all_cluster_freqs = h5['cluster_props'][:, 0][all_clusters]

    for i, feature in enumerate(features):
        label, mask = feature

        labels.append(label)
        counts[i] = np.sum(mask.astype(bool)[bins] * cluster_freqs)
        counts_all_clusters[i] = np.sum(
            mask.astype(bool)[all_bins] * all_cluster_freqs
        )

    plot_distribution(counts, labels, Y_comp=counts_all_clusters, figsize=figsize)


def plot_contact_logo(
    logo,
    feature,
    labels,
    colors,
    xtick_step_size=100,
    title: str = 'Contact Feature Logo',
    with_separate: bool = False,
):
    num_feature_states, num_bins = logo.shape

    extra_rows = num_feature_states if with_separate else 0

    fig, ax = plt.subplots(
        nrows=2 + extra_rows,
        ncols=1,
        figsize=(20, 6 + extra_rows),
        sharex=True,
        gridspec_kw=dict(
            height_ratios=([8] + [1] * extra_rows + [1]),
            hspace=0.01
        )
    )

    X = np.arange(num_bins)
    Y_min = logo.min()
    Y_max = logo.max()

    width = 1

    bottom = np.zeros(num_bins).astype(np.uint32)
    for state in np.arange(num_feature_states):
        ax[0].bar(
            X,
            logo[state],
            width,
            label=labels[state],
            bottom=bottom,
            color=colors[state]
        )
        bottom = bottom + logo[state]

    ax[0].set_title(title)
    ax[0].set_xticks([], [])
    ax[0].set_ylabel('Contact States')
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[0].legend()

    if extra_rows:
        for state in np.arange(num_feature_states):
            ax[state + 1].bar(
                X,
                logo[state],
                width,
                label=labels[state],
                color=colors[state]
            )
            ax[state + 1].set_ylim(Y_min, Y_max)
            ax[state + 1].set_xticks([], [])
            ax[state + 1].spines['top'].set_visible(False)
            ax[state + 1].spines['left'].set_visible(False)
            ax[state + 1].spines['right'].set_visible(False)
            ax[state + 1].spines['bottom'].set_visible(False)

    for state in np.arange(num_feature_states):
        ax[extra_rows + 1].bar(
            X,
            feature[state],
            width,
            label=labels[state],
            color=colors[state]
        )

    ax[extra_rows + 1].set_xticks(
        np.r_[np.arange(0, num_bins, xtick_step_size), num_bins - 1]
    )
    ax[extra_rows + 1].set_xlabel('Bins')
    ax[extra_rows + 1].set_yticks([],[])
    ax[extra_rows + 1].set_ylabel('State')
    ax[extra_rows + 1].spines['top'].set_visible(False)
    ax[extra_rows + 1].spines['left'].set_visible(False)
    ax[extra_rows + 1].spines['right'].set_visible(False)
    ax[extra_rows + 1].spines['bottom'].set_visible(False)
