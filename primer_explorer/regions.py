from primer_explorer.pcr import pcr_products_to_regions


class GenomeRegions:

    def __init__(self, bed_fhand=None, pcr_products=None):
        if bed_fhand is None and pcr_products is None:
            raise ValueError('At least one input is required')
        elif bed_fhand is not None and pcr_products is not None:
            raise ValueError('Can not work wth two inputs')
        self._bed_fhand = bed_fhand
        self._pcr_products = iter(pcr_products) if pcr_products is not None else None

    def __iter__(self):
        return self

    def __next__(self):
        if self._pcr_products is None:
            line_items = next(self._bed_fhand).split()
            chrom = line_items[0]
            start = int(line_items[1])
            stop = int(line_items[2])
        else:
            pcr_product = next(self._pcr_products)
            chrom = pcr_product[0].chrom_location[0]
            start = pcr_product[0].chrom_location[1]
            stop = pcr_product[1].chrom_location[1]
        return GenomeRegion(chrom, start, stop)


class GenomeRegion:

    def __init__(self, chrom, start, stop):
        self.chrom = chrom
        self.start = start
        self.stop = stop

    def __lt__(self, region2):
        if (self.chrom < region2.chrom or (self.chrom == region2.chrom and
                                           self.start < region2.start)):
            return True
        return False

        if self.chrom != region2.chrom:
            return True
        return self.start < region2.start

    def overlaps(self, region2):
        if self.chrom != region2.chrom:
            return False
        if self.stop <= region2.start:
            return False
        if self.start >= region2.stop:
            return False
        return True


def write_primer_regions_in_bed_format(product_results, out_dir, kmer_length,
                                       min_product_length=100,
                                       max_product_length=700,
                                       read_length=100):
    done_pairs = []
    for pcr_products_set in product_results:
        for pair, pcr_products in pcr_products_set['products'].items():
            if pair in done_pairs:
                continue
            out_path = out_dir / "{}-{}.bed".format(pair[0].decode(),
                                                    pair[1].decode())
            with out_path.open('w') as out_fhand:
                lines = pcr_products_to_regions(pcr_products, min_product_length,
                                                max_product_length, read_length,
                                                kmer_length=kmer_length)
                out_fhand.write("\n".join(lines) + "\n")
            done_pairs.append(pair)
