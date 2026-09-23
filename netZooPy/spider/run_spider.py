#!/usr/bin/env python

import sys
import getopt

from netZooPy.spider.spider import Spider


def main(argv):
    """
    Description:
        Run the SPIDER algorithm from the command line.

    Inputs:
        run_spider
        -h, --help       : help
        -e, --expression : expression values
        -m, --motif      : pair file of motif edges (TF, gene, weight)
        -f, --epifilter  : binary epigenetic filter matching the motif rows
        -p, --ppi        : pair file of PPI edges (TF, TF, weight)
        -o, --out        : output file
        -r, --rm_missing : remove missing (legacy mode only)

    Example:
        python run_spider.py -e expression.txt -m motif.txt \\
            -f epifilter.txt -p ppi.txt -o spider.txt

    Reference:
        Sonawane, Abhijeet Rajendra, et al. "Constructing gene regulatory
        networks using epigenetic data." npj Systems Biology and Applications
        7.1 (2021): 1-13.
    """
    expression_data = None
    motif = None
    epifilter = None
    ppi = None
    output_file = "output_spider.txt"
    rm_missing = False
    try:
        opts, args = getopt.getopt(
            argv,
            "he:m:f:p:o:r",
            ["help", "expression=", "motif=", "epifilter=", "ppi=", "out=", "rm_missing"],
        )
    except getopt.GetoptError:
        print(__doc__)
        sys.exit()
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(__doc__)
            sys.exit()
        elif opt in ("-e", "--expression"):
            expression_data = arg
        elif opt in ("-m", "--motif"):
            motif = arg
        elif opt in ("-f", "--epifilter"):
            epifilter = arg
        elif opt in ("-p", "--ppi"):
            ppi = arg
        elif opt in ("-o", "--out"):
            output_file = arg
        elif opt in ("-r", "--rm_missing"):
            rm_missing = arg

    if motif is None:
        print("Missing motif prior!")
        print(__doc__)
        sys.exit()

    print("Start SPIDER run ...")
    spider_obj = Spider(
        expression_data,
        motif,
        epifilter,
        ppi,
        save_tmp=True,
        remove_missing=rm_missing,
    )
    spider_obj.save_spider_results(output_file)
    print("All done!")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
