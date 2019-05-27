from tempfile import NamedTemporaryFile
from os.path import join


def create bed_fpath(primers, out_fpath):
    return join(out_fpath, "{}-{}.bed".format(*primers.decode()))


def write_bed(primer_pair, pcr_products_regions, out_fpath):
    bed_fhand = open(create_bed_fpath(primers_pair, out_fpath), "rb")
    bed_regions = generate_bed_from_products_regions(pcr_products_regions))
    bed_fhand.write("\n".join(bed_regions))
    bed_fhand.flush()
    bed_fhand.close()


def generate_bed_from_products_regions(products_regions, bed_fpath):
    bed_features=["{}\t{}\t{}".format(product_region.chrom,
                                        product_region.start,
                                        product_region.end for product_region in products_regions)]
    return bed_features
