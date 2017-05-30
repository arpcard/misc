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
            --max-target-seqs 1 --evalue 0.001 --threads 2 --max-hsps 1 "
    os.system('diamond blastp --in proteindb.fsa --db protein.db  \
            --query {} --outfmt 6 --out {}.diamond.txt \
            --salltitles {}'.format(i, o, cfilter))


def prepare(o):

    blast6table_header = "qseqid sseqid pident length mismatch gapopen \
                            qstart qend sstart send evalue bitscore".split()

    result_diamond = pd.read_csv('{}.diamond.txt'.format(o),
                                 header=None, names=blast6table_header,
                                 sep='\t')
    result_blast = pd.read_csv('{}.blast.txt'.format(o),
                               header=None, names=blast6table_header, sep='\t')

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

    # log1p transform to improve visualisation
    delta_log_evalue = np.log1p(b['evalue']) - np.log1p(d['evalue'])
    delta_bitscore = b['bitscore'] - d['bitscore']

    plot_graph(delta_log_evalue,
               'BLAST vs DIAMOND ln(E-value+1) Differences KDE',
               'diff_evalue_{}.png'.format(desc))

    plt.clf()
    sns.kdeplot(delta_log_evalue, delta_bitscore, kernel='gau')
    plt.title('BLAST vs DIAMOND bitscore and ln(E-value+1) Differences')
    plt.xlabel("delta(ln(E-value+1))")
    plt.ylabel("delta(bitscore)")
    plt.savefig('2d_evalue_bitscore_kde.png')

    plot_graph(delta_bitscore,
               'BLAST vs DIAMOND Bitscore Differences KDE',
               'diff_bitscore_{}.png'.format(desc))
    plot_graph(delta_bitscore / b['bitscore'],
               'BLAST vs DIAMOND Normalised Bitscore Differences KDE',
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
    if args.threaded:
        threaded(inputSeq, outputFile, desc)
    else:
        # non threaded
        non_threaded(inputSeq, outputFile, desc)

    print('Total running time {}s'.format(round(time.time() - t0, 3)))


def plot_graph(delta, title, filename='test.png'):
    plt.clf()
    sns.kdeplot(delta, kernel='gau')
    plt.title(title)
    plt.xlabel("delta scores")
    plt.ylabel("density")
    plt.savefig(filename)


def run():
    parser = argparse.ArgumentParser(description='benchmark tests')
    parser.add_argument('-i', '--input', dest="input", default=None,
                        required=True, help='input file (fasta)')
    parser.add_argument('-o', '--output', dest="output", default="summary",
                        required=True, help='output file (tab-delimited)')
    parser.add_argument('-t', '--threaded', dest="threaded",
                        action='store_false',
                        help='disable threading (default on)')
    parser.add_argument('-d', '--desc', dest="desc", default="na",
                        required=True, help='Description for the input')
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
