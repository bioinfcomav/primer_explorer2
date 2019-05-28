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


def filter_bed_by_annotation(bed_fpath, annotation_fpath):
    if bed_fpath.exists():
        bin = config.BEDTOOLS_INTERSECT_BINARY
        tmp_fhand = NamedTemporaryFile()
        args = ['-a', str(annotation_fpath),
                '-b', str(bed_fpath)]
        cmd = " ".join(bin + args).encode()
        results = _run_script(cmd, tmp_fhand, tmp_fhand)
        return results.decode('utf-8').rstrip()


def get_gff_stats(results):
    gff_stats = {}
    for library, genes in results.items():
        genes_with_exons = []
        mrnas = []
        five_prime_UTRs = []
        three_prime_UTRs = []
        for gene, features in genes.items():
            if len(features['exons']) > 0:
                genes_with_exons.append(gene)
            mrnas.append(gene)
            if "five_prime_UTR" in features:
                five_prime_UTRs.append(features["five_prime_UTR"])
            if "three_prime_UTR" in features:
                three_prime_UTRs.append(features["three_prime_UTR"])

        gff_stats[library] = {"mrnas_covered": len(mrnas),
                              "mrnas_with_exons_covered": len(genes_with_exons),
                              "genes_with_five_prime_utr": len(five_prime_UTRs),
                              "genes_with_three_prime_utr": len(three_prime_UTRs)}
    return gff_stats
