from collections import Counter

_DUST_CACHE = {}


def _rolling_window_serie(serie, window, length_, step):
    '''It yields lists of items with a window number of elements'''
    return (serie[i:i + window] for i in range(0, length_ - window + 1, step))


def _rolling_window_iter(iterator, window, step):
    '''It yields lists of items with a window number of elements giving
     an iterator'''
    items = []
    for item in iterator:
        if len(items) >= window:
            yield items
            items = items[step:]
        items.append(item)
    else:
        if len(items) >= window:
            yield items


def rolling_window(iterator, window, step=1):
    'It yields lists of items with a window number of elements'
    try:
        length_ = len(iterator)
    except TypeError:
        length_ = None
    if length_ is None:
        return _rolling_window_iter(iterator, window, step)
    else:
        return _rolling_window_serie(iterator, window, length_, step)


def _calculate_rawscore(string):
    'It returns a non-normalized dustscore'
    triplet_counts = Counter()
    for triplet in rolling_window(string, 3):
        # It should do something with non ATCG, but we sacrifice purity for
        # speed. Maybe we should reconsider this
        triplet_counts[triplet.upper()] += 1

    return sum(tc * (tc - 1) * 0.5 for tc in triplet_counts.values())


def calculate_dust_score(seq, windowsize, windowstep):
    '''It returns the dust score.
    From: "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA
    Sequences"
    doi:10.1089/cmb.2006.13.1028
    and re-implemented from PRINSEQ
    '''
    length = len(seq)
    if length == 3:
        return 0
    if length <= 5:
        return None

    dustscores = []
    if length > windowsize:
        windows = 0
        for seq_in_win in rolling_window(seq, windowsize, windowstep):
            score = _calculate_rawscore(seq_in_win)
            dustscores.append(score / (windowsize - 2))
            windows += 1
        remaining_seq = seq[windows * windowstep:]
    else:
        remaining_seq = seq

    if len(remaining_seq) > 5:
        length = len(remaining_seq)
        score = _calculate_rawscore(remaining_seq)
        dustscore = score / (length - 3) * (windowsize - 2) / (length - 2)
        dustscores.append(dustscore)

    # max score should be 100 not 31
    dustscore = sum(dustscores) / len(dustscores) * 100 / 31
    return dustscore


def dust_score_is_ok(seq, windowsize, windowstep, threshold):
    try:
        return _DUST_CACHE[seq]
    except KeyError:
        dust = calculate_dust_score(seq, windowsize, windowstep)
        result = dust < threshold
        _DUST_CACHE[seq] = result
        return result
