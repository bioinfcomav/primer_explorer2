import gzip


def get_fhand(fpath):
    if is_gzipfile(fpath):
        fhand = gzip.open(fpath, 'r')
    else:
        fhand = open(fpath, 'rb')
    return fhand


def is_gzipfile(fpath):
    with gzip.open(fpath) as fhand:
        try:
            fhand.readline()
            return True
        except IOError:
            return False
