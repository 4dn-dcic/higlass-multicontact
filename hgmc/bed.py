import gzip
import math
import numpy as np
import os
import sqlite3
import time

# bed_filepath needs to be sorted!
def bed_to_sql(bed_filepath, chrom_sizes_filepath, feature_name: str = None):
    base, ext = os.path.splitext(bed_filepath)
    open_bed = lambda f: open(f, 'r')

    if ext == '.gz' or ext == '.gzip':
        open_bed = lambda f: gzip.open(f, 'rt')
        base, ext = os.path.splitext(base)

    sqlite_filepath = base + ".sqlite"

    with open_bed(bed_filepath) as f:
        num_columns = len(f.readline().split('\t'))

        if num_columns < 4 and feature_name is None:
            raise ValueError(
                'Either provide a BED4 file or provide `feature_name`'
            )

        # Reset cursor
        f.seek(0)

        conn = sqlite3.connect(sqlite_filepath)

        conn.execute("DROP TABLE IF EXISTS intervals")
        conn.execute(
            """
            CREATE TABLE intervals
            (
                id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                chrom TEXT,
                start INT,
                end INT,
                name TEXT
            )
            """
        )

        conn.execute("DROP TABLE IF EXISTS features")
        conn.execute("CREATE TABLE features (name TEXT PRIMARY KEY NOT NULL)")

        cur = conn.cursor()

        def split_line_bed3(line):
            return (*line.split('\t'), feature_name)

        def split_line_bed4(line):
            return line.split('\t')[:4]

        split_line = split_line_bed4 if num_columns >= 4 else split_line_bed3

        unique_features = set()

        for line in f:
            chrom, start, end, name = split_line(line)

            unique_features.add(name)

            cur.execute(
                "INSERT INTO intervals (chrom, start, end, name) "
                "VALUES (?, ?, ?, ?)",
                (chrom, int(start), int(end), name),
            )

        conn.execute("CREATE INDEX interval_bounds ON intervals (start, end)")

        for name in unique_features:
            cur.execute("INSERT INTO features (name) VALUES (?)", (name,),)

        if chrom_sizes_filepath:
            conn.execute("DROP TABLE IF EXISTS chrom_sizes")

            conn.execute(
                """
                CREATE TABLE chrom_sizes
                (
                    chrom TEXT PRIMARY KEY,
                    size INT
                )
                """
            )

            cur = conn.cursor()

            with open(chrom_sizes_filepath, 'r') as f:
                for line in f:
                    chrom, size = line.split('\t')

                    cur.execute(
                        "INSERT INTO chrom_sizes (chrom, size) "
                        "VALUES (?, ?)",
                        (chrom, int(size)),
                    )

        conn.commit()
        conn.close()


def sql_features(sqlite_filepath):
    conn = sqlite3.connect(sqlite_filepath)
    cur = conn.cursor()

    return [x[0] for x in cur.execute("SELECT name FROM features").fetchall()]



def sql_coverage(
    sqlite_filepath,
    chrom: str,
    bin_size: int,
    features: list = None,
    # Number of bp of the feature that needs to be in the bin for the feature
    # to count
    count_at_feat_cov: int = None,
    # If true, `count_at_feat_cov` represents the percentage of the feature size
    # that needs to be in the bin for the feature to count
    rel_count_at_feat_cov: bool = False,
    # Number of bp that need to be cover by the feature for the feature
    # to count
    count_at_bin_cov: int = None,
     # If true, `count_at_bin_cov` represents the percentage of the bin size
     # that needs to be cover by the feature for the feature to count
    rel_count_at_bin_cov: bool = False,
    # By default, if `count_at_feat_cov` and `count_at_bin_cov` are specified
    # both conditions need to be fulfilled. If `feat_or_bin_cov` is true it's
    # enough that one is fulfilled
    feat_or_bin_cov: bool = False,
    timeit: bool = False,
    verbose: bool = False
):
    t0 = time.time()
    conn = sqlite3.connect(sqlite_filepath)
    cur = conn.cursor()

    res = cur.execute(
        "SELECT * FROM chrom_sizes WHERE chrom = ?", (chrom,),
    ).fetchone()

    if res is None:
        return None

    _, chrom_size = res

    res = cur.execute(
        "SELECT MIN(id), MAX(id) FROM intervals WHERE chrom = ?", (chrom,),
    ).fetchone()

    if res is None:
        return None

    min_id, max_id = res

    num_bins = math.ceil(chrom_size / bin_size)
    chrom_padded_size = bin_size * num_bins
    coverage = np.zeros(num_bins).astype(int)

    constraint_id = f'id >= {min_id} AND id <= {max_id} '
    constraint_feature = ''

    if features is not None:
        if isinstance(features, list):
            c = ' OR '.join([f'name = "{f}"' for f in features])
            constraint_feature = f'AND ({c}) '
        else:
            constraint_feature = f'AND name = "{features}" '

    i = 0
    if count_at_feat_cov is None and count_at_bin_cov is None:
        for bin_start in np.arange(0, chrom_padded_size, bin_size):
            bin_end = bin_start + bin_size
            count = cur.execute(
                f"""
                SELECT
                    COUNT(*)
                FROM
                    intervals
                WHERE
                    {constraint_id}
                    {constraint_feature}
                    AND start < ?
                    AND end >= ?
                """,
                (int(bin_end), int(bin_start),),
            ).fetchone()[0]

            coverage[i] = count
            if verbose:
                print(f'In [{bin_start},{bin_end}) found {res[0]} annotations')

            i += 1
    else:
        for bin_start in np.arange(0, chrom_padded_size, bin_size):
            bin_end = bin_start + bin_size
            results = cur.execute(
                f"""
                SELECT
                    start, end, end-start
                FROM
                    intervals
                WHERE
                    {constraint_id}
                    {constraint_feature}
                    AND start < ?
                    AND end >= ?
                """,
                (int(bin_end), int(bin_start),),
            ).fetchall()

            if results is not None:
                count = 0
                for result in results:
                    feat_start, feat_end, feat_size = result
                    feat_coverage = min(feat_end, bin_end) - max(feat_start, bin_start)

                    should_count = True
                    if count_at_feat_cov:
                        threshold = rel_count_at_feat_cov * count_at_feat_cov * feat_size or count_at_feat_cov
                        should_count = feat_coverage >= threshold
                        if feat_or_bin_cov:
                            if should_count:
                                count += 1
                                continue
                            else:
                                should_count = True

                    if should_count and count_at_bin_cov:
                        threshold = rel_count_at_bin_cov * count_at_bin_cov * bin_size or count_at_bin_cov
                        should_count = feat_coverage >= threshold

                    if should_count:
                        count += 1

                coverage[i] = count
                if verbose:
                    print(f'In [{bin_start},{bin_end}) found {res[0]} annotations')

            i += 1

    if timeit:
        print(f'Took {(time.time() - t0):.3f} sec')

    return coverage
