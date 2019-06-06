import subprocess
from collections import Counter
from tempfile import NamedTemporaryFile

from primer_explorer import config


def count_features_in_gff(gff_fhand, features_to_count=None):

    feature_counts = Counter()
    for line in gff_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        line = line.split("\t")

        feature = line[2]
        if (features_to_count is None or
                (features_to_count and feature in features_to_count)):
            feature_counts[feature] += 1
    return feature_counts


def count_regions_in_gff_features(bed_fpath, gff_fpath,
                                  features_to_count=None):
    bedtools_binary = config.BEDTOOLS_BINARY
    cmd = [bedtools_binary, 'intersect', "-a", str(gff_fpath), "-b",
           str(bed_fpath)]

    result_fhand = NamedTemporaryFile('wb')
    subprocess.run(cmd, stdout=result_fhand, check=True)
    result_fhand.flush()
    with open(result_fhand.name) as gff_fhand:
        counts = count_features_in_gff(gff_fhand,
                                       features_to_count=features_to_count)
    result_fhand.close()
    return counts
