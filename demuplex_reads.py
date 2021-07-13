#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
author: Xiaoping
date: 06/27/2019


take primer sequences and match to the beginning of the sequence, if the mismtach <= n, then that sequence can be considered into that primer


modified: 07/24/2019
modified-content: modified that --reads parameter can take in more reads

modified: 07/31/2019
modified-content: added mismatch number to matching primers
"""

import os
import gzip
import re
import argparse
import sys
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

__author__ = "Xiaoping Li"
__version__ = "0.0.1"
__email__ = "lixiaopi@oregonstate.edu"





def getInput():

    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", dest="readpath", nargs="+", help="raw read path")
    parser.add_argument("--m", dest="mismatch", help="mismatch number")
    parser.add_argument("--primer", dest="primer", help="primer file")


    args = parser.parse_args()
    readp = args.readpath # alist
    mismatch = args.mismatch
    primer = args.primer

    return readp, primer, mismatch


# def subFunc(x, *args):
#     return re.sub(args[0], args[1], x)


def mismatchStr(str1, str2):
    degen_primers = {"Y":["C", "T"],
                    "M":["A", "C"],
                    "W":["A", "T"],
                    "R":["A", "G"],
                    "S":["C", "G"],
                    "K":["G", "T"]}


    count = 0
    for i, v in enumerate(str1):

        if v in degen_primers.keys():
            if not str2[i] in degen_primers[v]:
                count += 1
        else:
            if str2[i] != v:
                count += 1

    return count



def split_reads(parser):
    """[summary]
    
    Arguments:
        parser {tuple} -- from getInput()
        n {int} -- mismatch
    
    Raises:
        FileNotFoundError: [description]
        FileNotFoundError: [description]
    
    Returns:
        None -- Split the reads into category defined by primerfile
    """

    raw_reads_path = parser[0] # a list
    primer_file = parser[1] # a list
    n = parser[2]


    if not os.path.exists(primer_file):
        raise FileNotFoundError(f"Can not find primer file {raw_reads_path}")


    primers = pd.read_csv(primer_file, delimiter="\t")

    subfolders = primers.folder

    primerF = primers.primer_f # get only forward primer

    # primerF = primers.primer_f.apply(subFunc, args=(pattern_degen, "."))
    # primerR = primers.primer_r.apply(subFunc, args=(pattern_degen, "."))

    for i in subfolders:
        if not os.path.exists(i):
            os.makedirs(i)

    folder_abspath = [os.path.abspath(i) for i in subfolders]

    # create mapping for folder and primer pairs: {"ITSpath":"^f|^r"}

    if raw_reads_path:
        list_reads = [ os.path.join(i, j) for i in raw_reads_path for j in os.listdir(i)]

        print(f"Total reads {len(list_reads)}")

    else:
        raise FileNotFoundError(f"Can not find raw reads collection{raw_reads_path}")


    # get forward reads and reverse reads
    pattern_r1 = re.compile('.*_R1_.*')
    pattern_r2 = re.compile('.*_R2_.*')

    R1 = [ i for i in list_reads if re.match(pattern_r1, os.path.basename(i))]
    R2 = [ i for i in list_reads if re.match(pattern_r2, os.path.basename(i))]
    R1.sort(key=lambda x: int(os.path.basename(x).split("-")[1].split("s")[1]))
    R2.sort(key=lambda x: int(os.path.basename(x).split("-")[1].split("s")[1]))
    R1R2_list = list(zip(R1, R2))



    for index, val in enumerate(folder_abspath):

        primer_selected = primerF[index]

        for reads in tqdm(R1R2_list, desc=subfolders[index]):
            f_read = reads[0]
            r_read = reads[1]

            with gzip.open(f_read, 'rt') as handle_f, gzip.open(r_read, 'rt') as handle_r, gzip.open(os.path.join(val, os.path.basename(f_read)), 'wt') as output_f, gzip.open(os.path.join(val, os.path.basename(r_read)), 'wt') as output_r:

                # get forward reads
                record_f = list(SeqIO.parse(handle_f, 'fastq'))
                record_r = list(SeqIO.parse(handle_r, 'fastq'))

                for i, v in enumerate(record_f):

                    seq = str(v.seq)
                    
                    # match string
                    mismatch = mismatchStr(primer_selected, seq[0:len(primer_selected)])

                    if mismatch <= int(n):
                        SeqIO.write(v, output_f, 'fastq')
                        SeqIO.write(record_r[i], output_r, 'fastq')

            
    print("Completed!")
    return "Completed"
    




if __name__ == "__main__":

    parser = getInput()
    split_reads(parser)

