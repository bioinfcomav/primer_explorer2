import os


def create_bed_fpath(primers, out_fpath):
    return os.path.join(out_fpath, "{}-{}.bed".format(primers[0].decode(),
                                                      primers[1].decode()))


def write_bed(primers_pair, pcr_products_regions, out_fpath):
    bed_fhand = open(create_bed_fpath(primers_pair, out_fpath), "wb")
    bed_regions = generate_bed_from_products_regions(pcr_products_regions)
    bed_fhand.write("\n".join(bed_regions).encode())
    bed_fhand.flush()
    bed_fhand.close()


def generate_bed_from_products_regions(products_regions):
    bed_features = []
    for product_region in products_regions:
        bed_features.append("{}\t{}\t{}".format(product_region.chrom.decode(),
                                                product_region.start,
                                                product_region.stop))
    return bed_features
