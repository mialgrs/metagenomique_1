#! /bin/env python3
# -*- coding: utf-8 -*-

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import textwrap
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Mia Legras"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Mia Legras"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mia Legras"
__email__ = "mialegras@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file',
        type=isfile, required=True,
        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen',
        type=int, default = 400,
        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount',
        type=int, default = 10,
        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size',
                        type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size',
                        type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen): # attention : yield marche pas avec fasta normaux
    isfile(amplicon_file)
    with gzip.open(amplicon_file, 'rt') as fasta:
        seq = ""
        for line in fasta:
            if line.startswith('>'):
                if len(seq) > minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line.strip() #recup ligne suiv
        if len(seq) >= minseqlen:
            yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    seq = list(read_fasta(amplicon_file, minseqlen))
    uniqseq = set(seq)
    res = []
    for subseq in uniqseq:
        occur = seq.count(subseq)
        if occur >= mincount:
            res.append([subseq, occur])
    res = sorted(res, key=lambda x : x[1], reverse=True)
    for i in res:
        yield i


def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    seq1 = list(alignment_list[0])
    seq2 = list(alignment_list[1])
    nb_nucleo_id = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            nb_nucleo_id += 1
    identity = (nb_nucleo_id / len(seq1))*100
    return identity


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount,
                                chunk_size, kmer_size):
    seq = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    OTU = [seq[0]]
    for seq1 in seq[1:]:
        for seq2 in OTU:
            align = nw.global_align(seq1[0], seq2[0], gap_open=-1,
                gap_extend=-1, matrix=os.path.abspath(os.path.join(
                os.path.dirname(__file__),"MATCH")))
        if get_identity(align)<97:
            OTU.append(seq1)
            print(OTU)
    return OTU


def write_OTU(OTU_list, output_file):
    with open(output_file, 'w') as file:
        for i in range(len(OTU_list)):
            file.write(f'>OTU_{i+1} occurrence:{OTU_list[i][1]}\n')
            file.write(f'{textwrap.fill(OTU_list[i][0], width=80)}\n')

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, 
        args.mincount)
    print(otu)


if __name__ == '__main__':
    main()
