from collections import Counter, defaultdict, OrderedDict
from os import listdir
from os.path import join
from subprocess import CalledProcessError, check_output
from tempfile import NamedTemporaryFile


from primer_explorer import config


def _run_script(cmd):
    tmp_fhand = NamedTemporaryFile()
    tmp_fhand.write(cmd)
    tmp_fhand.flush()
    try:
        return check_output(['sh', tmp_fhand.name])
    except CalledProcessError:
        raise


def count_bedtools_intersect_features(stdout):
    features = Counter()
    stdout = stdout.decode('utf-8').rstrip()
    if len(stdout) > 0:
        for line in stdout.split("\n"):
            line = line.split("\t")
            feature = line[2]
            if feature in config.GFF_FEATURES:
                features[feature] += 1
    return features


def get_gff_intersect_results(bed_fpaths, gff, gff_results):
    for primer_pairs, bed_fpath in bed_fpaths.items():
        bin = config.BEDTOOLS_INTERSECT_BINARY
        args = ["-a", str(gff), "-b", str(bed_fpath)]
        cmd = " ".join(bin + args).encode()
        stdout = _run_script(cmd)
        results = count_bedtools_intersect_features(stdout)
        gff_results[primer_pairs].update(results)
    return gff_results
