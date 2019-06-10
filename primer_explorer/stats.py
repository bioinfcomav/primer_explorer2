from collections import Counter


from primer_explorer.primer3.primer3 import reverse_complement

LABELS = {'title': 'histogram', 'xlabel': 'values',
          'ylabel': 'count', 'minimum': 'minimum',
          'maximum': 'maximum', 'average': 'average',
          'variance': 'variance', 'sum': 'sum',
          'items': 'items', 'quartiles': 'quartiles'}

PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES = "percentage_sequenciable_nucleotides"
ADJUSTED_PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES = "adjusted_percentage_sequenciable_nucleotides"
NUM_OF_POSSIBLE_PRODUCTS_10000 = "num_possible_products_10000"
NUM_OF_POSSIBLE_PRODUCTS_700 = "num_possible_products_700"
NUM_SEQUENCIABLE_PRODUCTS = "num_sequenciable_products"
NUM_SHORT_PRODUCTS = "num_short_products"
NUM_UNION_SITES = "num_union_sites"
NUM_REPETITIVE_REPETITIVE_PRODUCTS = "num_repetitive-repetitive_products"
NUM_UNIQUE_UNIQUE_PRODUCTS = "num_unique-unique_products"
NUM_UNIQUE_REPETITIVE_PRODUCTS = "num_unique-repetitive_products"


class IntCounter(Counter):
    '''This is a subclass  of the Counter on python collections.

    It adds some statistical functionalities to the class'''

    def __init__(self, *args, **kwargs):
        'It inialices a '
        super(IntCounter, self).__init__(*args, **kwargs)
        self.labels = LABELS.copy()

    @property
    def min(self):
        'Get the minimun value'
        return min(self.keys())

    @property
    def max(self):
        'Get the maximun value'
        return max(self.keys())

    @property
    def count(self):
        'It returns the count of the values stored in the array'
        return sum(self.values())

    @property
    def median(self):
        'It calculates the median of the values appended'
        quotient, remainder = divmod(self.count, 2)
        if remainder == 0:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            return (val1 + val2) / 2
        else:
            return self._get_value_for_index(quotient)

    @property
    def sum(self):
        'It gets the sum of the values'
        sum_ = 0
        for index, value in self.items():
            sum_ += (index * value)
        return int(sum_)

    @property
    def average(self):
        'It calculates the average'
        count = self.count
        sum_ = self.sum
        return sum_ / count

    @property
    def variance(self):
        'It gets the variance of the values'
        mean = self.average
        sum_ = 0
        for index, counts in self.items():
            sum_ += ((index - mean) ** 2) * counts
        return sum_ / self.count

    @property
    def quartiles(self):
        'It returns the quartiles'
        num_items = self.count
        if num_items < 4:
            msg = 'At least 4 values are required to calculate the quartiles'
            raise RuntimeError(msg)
        # quartile 1
        quotient, remainder = divmod(num_items + 1, 4)
        if not remainder:
            quartile1 = self._get_value_for_index(quotient - 1)
        else:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            quartile1 = (val1 + val2) / 2
        # quartile 3
        quotient, remainder = divmod((num_items + 1) * 3, 4)
        if not remainder:
            quartile3 = self._get_value_for_index(quotient - 1)
        else:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            quartile3 = (val1 + val2) / 2
        return quartile1, self.median, quartile3

    @property
    def irq(self):
        'It gets the interquartile range'
        # pylint: disable=W0612
        quart1 = self.quartiles[0]
        quart3 = self.quartiles[2]
        return quart3 - quart1

    @property
    def outlier_limits(self):
        'It returns the intercuartile'
        # pylint: disable=W0612
        quart1 = self.quartiles[0]
        quart3 = self.quartiles[2]

        iqr = self.irq
        limit_distance = round(iqr * 1.5)

        start = int(quart1 - limit_distance)
        end = int(quart3 + limit_distance)
        return (start, end)

    def _get_value_for_index(self, position):
        '''It takes a position and it returns the value for the given index'''
        cum_count = 0
        for index in sorted(self.keys()):
            count = self[index]
            cum_count += count
            if position <= cum_count - 1:
                return index
        else:
            if position >= cum_count:
                raise IndexError('You asked for an index beyond the scope')
            return index

    def _calculate_dist_range(self, min_, max_, outlier_threshold):
        'it calculates the range for the histogram'
        if ((min_ is not None or max_ is not None) and
                outlier_threshold is not None):
            msg = 'You can not pass max, min and outlier_threslhosld to '
            msg += 'calculate distribution range'
            raise ValueError(msg)

        if min_ is None:
            min_ = self.min
        if max_ is None:
            max_ = self.max
        if outlier_threshold:
            left_limit = self.count * outlier_threshold / 100
            rigth_limit = self.count - left_limit
            left_value = self._get_value_for_index(left_limit)
            rigth_value = self._get_value_for_index(rigth_limit)

            if min_ < left_value:
                min_ = left_value
            if max_ > rigth_value:
                max_ = rigth_value
        return min_, max_

    def calculate_bin_edges(self, min_, max_, n_bins=None):
        'It calculates the bin_edges'
        min_bins = 10
        max_bins = 20
        if n_bins is None:
            num_values = int(max_ - min_)
            if num_values == 0:
                n_bins = 1
            elif num_values < min_bins:
                n_bins = num_values
            else:
                n_bins = int(self.count / 10)
                if n_bins < min_bins:
                    n_bins = min_bins
                if n_bins > max_bins:
                    n_bins = max_bins
                if n_bins > num_values:
                    n_bins = num_values

        # now we can calculate the bin edges
        distrib_span = max_ - min_ if max_ != min_ else 1

        if distrib_span % n_bins:
            distrib_span = distrib_span + n_bins - (distrib_span % n_bins)
        bin_span = distrib_span // n_bins
        bin_edges = [min_ + bin_ * bin_span for bin_ in range(n_bins + 1)]
        return bin_edges

    def calculate_distribution(self, bins=None, min_=None, max_=None,
                               outlier_threshold=None):
        'It returns an histogram with the given range and bin'
        if self.count == 0:
            raise RuntimeError('No items in IntCounter')
        distrib = []
        min_, max_ = self._calculate_dist_range(min_, max_, outlier_threshold)
        if min_ is None or max_ is None:
            return None
        bin_edges = self.calculate_bin_edges(min_, max_, bins)
        for bin_index, left_edge in enumerate(bin_edges):
            try:
                rigth_edge = bin_edges[bin_index + 1]
            except IndexError:
                break
            sum_values = 0

            for index2 in sorted(self.keys()):
                value = self[index2]
                if index2 > rigth_edge:
                    break

                elif (left_edge <= index2 and index2 < rigth_edge or
                      left_edge <= index2 and index2 == max_):
                    sum_values += value

            distrib.append(sum_values)
        return {'counts': distrib, 'bin_limits': bin_edges}

    def update_labels(self, labels):
        'It prepares the labels for output files'
        self.labels.update(labels)

    def count_relative_to_value(self, value, comparison):
        'It counts the ints greater, equal, etc, relative to the given value.'
        return sum([counts for val, counts in self.items()
                    if comparison(val, value)])

    def __add__(self, other):
        'Add counts from two counters.'
        counter_python = super(IntCounter, self).__add__(other)
        return self.__class__(counter_python)

    def __str__(self):
        'It writes some basic stats of the values'
        if self.count != 0:
            labels = self.labels

            # now we write some basic stats
            def format_num(x):
                return '{:,d}'.format(x) if isinstance(x, int) else '%.2f' % x

            text = '{}: {}\n'.format(labels['minimum'], format_num(self.min))
            text += '{}: {}\n'.format(labels['maximum'], format_num(self.max))
            text += '{}: {}\n'.format(labels['average'],
                                      format_num(self.average))

            if labels['variance'] is not None:
                text += '{}: {}\n'.format(labels['variance'],
                                          format_num(self.variance))
            if labels['sum'] is not None:
                text += '{}: {}\n'.format(labels['sum'],
                                          format_num(self.sum))
            if labels['items'] is not None:
                text += '{}: {}\n'.format(labels['items'], self.count)
            text += '\n'
            return text
        return ''


def get_total_nondimer_pcr_products(pcr_products):
    non_dimer_pcr_products = []
    for product in pcr_products:
        fwd_primer = product[0].seq
        rev_primer = product[1].seq
        if fwd_primer == reverse_complement(rev_primer):
            continue
        else:
            non_dimer_pcr_products.append(product)
    return non_dimer_pcr_products


def filter_by_length(pcr_products, min_length=None, max_length=None):
    viable_products = []
    for product in pcr_products:
        fwd_primer_position = product[0].chrom_location[1]
        rev_primer_position = product[1].chrom_location[1]
        length = abs(fwd_primer_position - rev_primer_position)
        min_ok = False
        if (min_length is None or length >= min_length):
            min_ok = True
        max_ok = False
        if (max_length is None or length <= max_length):
            max_ok = True
        if min_ok and max_ok:
            viable_products.append(product)

    return viable_products


def get_union_sites_for_primers(pcr_products):
    union_sites = Counter()
    for product in pcr_products:
        fwd_primer = product[0].seq
        rev_primer = product[1].seq
        union_sites[fwd_primer] += 1
        union_sites[rev_primer] += 1
    return union_sites


def get_total_euchromatic_products(pcr_products):
    euchromatic_products = []
    for product in pcr_products:
        if all(not primer.is_heterochromatic for primer in product):
            euchromatic_products.append(product)
    return euchromatic_products


def get_total_heterochromatic_products(pcr_products):
    heterochromatic_products = []
    for product in pcr_products:
        if all(primer.is_heterochromatic for primer in product):
            heterochromatic_products.append(product)
    return heterochromatic_products


def get_total_mixed_products(pcr_products):
    mixed_products = []
    for product in pcr_products:
        if (not all(not primer.is_heterochromatic for primer in product) and
                not all(primer.is_heterochromatic for primer in product)):
            mixed_products.append(product)
    return mixed_products


def get_pcr_products_counts(pcr_products, min_length, max_length, genome_length):
    pcr_products_counts = Counter()
    total_products = get_total_nondimer_pcr_products(pcr_products)
    pcr_products_counts[NUM_OF_POSSIBLE_PRODUCTS_10000] = len(total_products)
    union_sites = get_union_sites_for_primers(pcr_products)
    pcr_products_counts[NUM_UNION_SITES] = union_sites

    sequenciable_products = filter_by_length(total_products, min_length, max_length)
    pcr_products_counts[NUM_SEQUENCIABLE_PRODUCTS] = {'min': min_length,
                                                      'max': max_length,
                                                      'count': len(sequenciable_products)}
    short_product_count = filter_by_length(total_products, max_length=min_length)
    pcr_products_counts[NUM_SHORT_PRODUCTS] = {'max': min_length,
                                               'count': len(short_product_count)}
    pcr_products_counts[NUM_OF_POSSIBLE_PRODUCTS_700] = len(sequenciable_products) + len(short_product_count)
    euchromatin_products = get_total_euchromatic_products(sequenciable_products)
    pcr_products_counts[NUM_UNIQUE_UNIQUE_PRODUCTS] = len(euchromatin_products)
    heterochromatin_products = get_total_heterochromatic_products(sequenciable_products)
    pcr_products_counts[NUM_REPETITIVE_REPETITIVE_PRODUCTS] = len(heterochromatin_products)
    mixed_products = get_total_mixed_products(sequenciable_products)
    pcr_products_counts[NUM_UNIQUE_REPETITIVE_PRODUCTS] = len(mixed_products)
    # TODO JUNTO BREADTH
    pcr_products_counts[PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES] = get_pcr_nucleotide_count(sequenciable_products,
                                                                                           genome_length)
    pcr_products_counts[ADJUSTED_PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES] = get_pcr_nucleotide_count(sequenciable_products,
                                                                                                    genome_length,
                                                                                                    scale_to_illumina_sequences=True)
    return pcr_products_counts


def get_stats_by_pair_in_sets(products_sets, min_length=100, max_length=1000):
    stats = {}
    for set_index, set_info in enumerate(products_sets):
        stats[set_index] = {'primers': set_info['primers'],
                            'stats': {}}
        for pair in sorted(set_info['products'].keys()):
            pcr_products = set_info['products'][pair]
            counts = get_pcr_products_counts(pcr_products,
                                             min_length=min_length,
                                             max_length=max_length)
            stats[set_index]['stats'][pair] = counts
    return stats


def get_pcr_nucleotide_count(pcr_products, genome_length,
                             scale_to_illumina_sequences=False,
                             illumina_min_pair_length=300):
    nucleotides = 0
    for product in pcr_products:
        start = product[0].chrom_location[1]
        end = product[1].chrom_location[1]
        distance = abs(start - end)
        if scale_to_illumina_sequences and distance > illumina_min_pair_length:
            distance = illumina_min_pair_length
        nucleotides += distance
    return float(nucleotides / genome_length)
