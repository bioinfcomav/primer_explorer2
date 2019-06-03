import subprocess
from collections import Counter

from primer_explorer import config


def count_features_in_gff(gff_content, features_to_count=None):
    feature_counts = Counter()
    for line in gff_content:
        line = line.strip()
        if not line or line.startwith('#'):
            continue
        line = line.split("\t")
        feature = line[2]
        if (features_to_count is None or
                (features_to_count and feature in config.GFF_FEATURES)):
            feature_counts[feature] += 1
    return feature_counts


def count_regions_in_gff_features(bed_fpath, gff_fpath,
                                  features_to_count=None):
    bedtools_binary = config.BEDTOOLS_BINARY
    cmd = [bedtools_binary, 'intersect', "-a", str(gff_fpath), "-b",
           str(bed_fpath)]
    output = subprocess.run(cmd, capture_output=True)
    if output:
        counts = count_features_in_gff(output,
                                       features_to_count=features_to_count)
    else:
        counts = None
    return counts
