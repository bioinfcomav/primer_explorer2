class GenomeRegions:

    def __init__(self, bed_fhand):
        self._bed_fhand = bed_fhand

    def __iter__(self):
        return self

    def __next__(self):
        line_items = next(self._bed_fhand).split()
        chrom = line_items[0]
        start = int(line_items[1])
        stop = int(line_items[2])
        return GenomeRegion(chrom, start, stop)


class GenomeRegion:

    def __init__(self, chrom, start, stop):
        self.chrom = chrom
        self.start = start
        self.stop = stop

    def __lt__(self, region2):
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