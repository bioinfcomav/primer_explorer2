#!/usr/bin/env python

import pickle
import argparse

from primer_explorer.stats import write_report


def parse_arguments():
    desc = "write stats about the searched primer combinations"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-p', '--pcr_products', required=True,
                        help='Path to write pcr_products result',
                        type=argparse.FileType('rb'))
    parser.add_argument('-o', '--output', required=True,
                        help='Path to write the report',
                        type=argparse.FileType('wb'))
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    products_fhand = args.pcr_products
    report_fhand = args.report_fhand

    return {'products_fhand': products_fhand, 'report_fhand': report_fhand}


def main():
    args = get_args()
    pcr_products_sets = pickle.load(args['products_fhand'])
    write_report(args['report_fhand'], pcr_products_sets)


if __name__ == "__main__":
    main()
