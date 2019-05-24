def generate_pair_stats(primer_pair, count):
    report = []
    report.append("{},{}".format(primer_pair[0].decode(), primer_pair[1].decode()))
    report.append("-" * 20)
    report.append("Total Number of products:\t{}".format(str(count["total_products"])))
    report.append("Number of viable pcr products:\t{}".format(str(count["viable_products"])))
    report.append("Number of nondimer and viable products (effective products):\t{}".format(str(count["total_non_dimer_products"])))
    report.append("Number of euchromatic effective products:\t{}".format(str(count["euchromatin_products"])))
    report.append("Number of heterochromatic effective products:\t{}".format(str(count["heterochromatin_products"])))
    report.append("Number of mixed effective products:\t{}".format(str(count["mixed_products"])))
    report.append("Ratio euchromatin (without mixed):\t{:.2f}".format(count["euchromatin_products"] / (count["viable_products"] - count["mixed_products"])))
    report.append("Ratio heterochromatin (without mixed):\t{:.2f}".format(count["heterochromatin_products"] / (count["viable_products"] - count["mixed_products"])))
    return report


def write_detailed_report(report_fhand, stats):
    report = []
    for idx, set_stats in stats.items():
        primers = set_stats['primers']
        report.append("PRIMER SET {}".format(str(idx)))
        report.append('Primers: ' + ', '.join([p.decode() for p in primers]))
        report.append("#" * 30)
        for pair, counts in set_stats['stats'].items():
            report.extend(generate_pair_stats(pair, counts))
            report.append("-" * 20)

    report_fhand.write("\n".join(report))
    report_fhand.flush()
