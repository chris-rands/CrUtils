from __future__ import division

import sys
import argparse

'''
*"N50 statistic is a measure of the average length of a set of sequences, with greater weight given to longer sequences.
It is used widely in genome assembly, especially in reference to contig lengths within a draft assembly. Given a set of sequences of varying lengths,
the N50 length is defined as the length N for which half of all bases in the sequences are in a sequence of length L < N.
In other words, N50 is the contig length such that using equal or longer contigs produces half the bases of the genome."--Wikipedia
'''

def get_args():
    parser = argparse.ArgumentParser(description='Calculate the Calculate N50 statistic* for set of contigs/scaffolds/chromosomes @Chris Rands 2016')
    parser.add_argument('-i', '--input',required=True, help='input fasta file name, use "-" for stdin')
    return parser.parse_args()

def yield_lines(fasta_file):
    if fasta_file == '-':
        for line in sys.stdin:
            yield line
    else:
        with open(fasta_file) as f:
            for line in f:
                yield line
    
def calculate_N50(gen):
    scaffold_sizes, scaffold_len = [], None
    for line in gen:
        if line == "\n":
            continue
        if line.startswith(">"):
            if scaffold_len is not None:
                scaffold_sizes.append(scaffold_len)
            scaffold_len = 0
            continue
        scaffold_len += len(line.rstrip())
    scaffold_sizes.append(scaffold_len)
    counter = 0
    threshold = sum(scaffold_sizes)/2
    for scaffold_size in sorted(scaffold_sizes):
        if (counter + scaffold_size) >= threshold:
            return scaffold_size
        else:
            counter += scaffold_size
    assert False, 'How did you get here?'

if __name__ == '__main__':
    N50 = calculate_N50(yield_lines(get_args().input))
    print("N50 {}bp\nN50 {}kb".format(N50,N50/1000))
