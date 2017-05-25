import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import csv
import multiprocessing
import time


def worker(name, i, o):
    """thread worker function"""
    if name == "diamond":
        diamond(i, o)
    else:
        blast(i, o)


def blast(i, o):
    os.system('blastp -query {} -db protein.db -evalue 0.001 \
              -max_target_seqs 1 -max_hsps 1 -num_threads 32 \
              -outfmt 6 -out {}.blast.txt'.format(i, o))


def diamond(i, o):
    cfilter = "--index-chunks 1 --block-size 1 --quiet --more-sensitive \
            --max-target-seqs 1 --evalue 0.001 --threads 32 --max-hsps 1 "
    os.system('diamond blastp --in proteindb.fsa --db protein.db  \
            --query {} --outfmt 6 --out {}.diamond.txt \
            --salltitles {}'.format(i, o, cfilter))


def prepare(o):
    result_diamond = pd.read_csv('{}.diamond.txt'.format(o))
    result_blast = pd.read_csv('{}.blast.txt'.format(o))

    return (result_blast, result_diamond)


def async_blast_diamond(i, o):
    p1 = multiprocessing.Process(target=worker, args=("blast", i, o))
    p2 = multiprocessing.Process(target=worker, args=("diamond", i, o))
    p1.start()
    p2.start()


def async_prepare_plot(o, desc):
    b, d = prepare(o)
    plot(b, d, desc)


def threaded(i, o, desc):
    """ Run one thread for blast then one thread for prepare and plot """
    blast_thread = multiprocessing.Process(target=async_blast_diamond,
                                           args=(i, o))
    blast_thread.start()
    blast_thread.join()

    prepare_plot_thread = multiprocessing.Process(target=async_prepare_plot,
                                                  args=(o, desc))
    prepare_plot_thread.start()


def non_threaded(i, o, desc):
    """ blast and diamond """
    # diamond
    diamond(i, o)
    # blast
    blast(i, o)

    """ prepare and plot """
    b, d = prepare(o)
    plot(b, d, desc)


def plot(b, d, desc):
    plot_graph(b['evalue'] - d['evalue'],
               'diff_evalue_{}.png'.format(desc))
    plot_graph((b['evalue'] - d['evalue']) / b['evalue'],
               'norm_evalue_{}.png'.format(desc))

    plot_graph(b['bitscore'] - d['bitscore'],
               'diff_bitscore_{}.png'.format(desc))
    plot_graph((b['bitscore'] - d['bitscore']) / b['bitscore'],
               'norm_bitscore_{}.png'.format(desc))


def main(args):
    t0 = time.time()
    if args.input is None:
        exit("Error - Missing input file")

    # inputs
    inputSeq = args.input
    outputFile = args.output
    desc = args.desc

    # threaded
    if args.threaded.lower() == "yes":
        threaded(inputSeq, outputFile, desc)
    else:
        # non threaded
        non_threaded(inputSeq, outputFile, desc)

    print('Total running time {}s'.format(round(time.time() - t0, 3)))


def plot_graph(delta, filename='test.png'):
    plt.clf()
    sns.kdeplot(delta, kernel='gau')
    plt.title("Gaussian Kernel Density Estimation \
            - DIAMOND vs BLAST bitscores")
    plt.xlabel("delta scores")
    plt.ylabel("density")
    plt.savefig(filename)


def run():
    parser = argparse.ArgumentParser(description='benchmark tests')
    parser.add_argument('-i', '--input', dest="input", default=None,
                        help='input file (fasta)')
    parser.add_argument('-o', '--output', dest="output", default="summary",
                        help='output file (tab-delimited)')
    parser.add_argument('-t', '--threaded', dest="threaded", default="yes",
                        help='enabled threading. Options are Yes or No \
                                (Default: Yes)')
    parser.add_argument('-d', '--desc', dest="desc", default="na",
                        help='Description for the input')
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
