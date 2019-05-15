from primer_explorer import config
from collections import Counter

_GC_RESULT_CACHE = {}


def blacklisted_seqs_in_seq(kmer):
    if b'N' in kmer:
        return True
    if b'AAAA' in kmer:
        return True
    if b'GGGG' in kmer:
        return True
    if b'TTTT' in kmer:
        return True
    if b'CCCC' in kmer:
        return True


def calculate_gc(seq):
    count = Counter(seq)
    g_count = count[config.G_CHAR]
    c_count = count[config.C_CHAR]
    gc = (g_count + c_count) / len(seq)
    return gc


def gc_is_ok(kmer, min_gc, max_gc):
    try:
        return _GC_RESULT_CACHE[kmer]
    except KeyError:
        if min_gc <= calculate_gc(kmer) <= max_gc:
            result = True
        else:
            result = False
        _GC_RESULT_CACHE[kmer] = result
        return result
