#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
author: Xiaoping
date: 06/27/2019
"""
import os
import pandas as pd
import re
from collections import defaultdict
import shutil
from tqdm import tqdm

class Sort_reads:

    def __init__(self, readpath, groupfile):
        """[summary]
        
        Arguments:
            readpath {[string]} -- reads path
            groupfile {[string]} -- usually your metafile (csv)
            outputfolder {[string]} -- the base directory to store sorted reads
        """
        self.readpath = readpath
        self.groupfile = groupfile
        self.grouped_reads = defaultdict(list)

    def __repr__(self):
        return "sort_reads()"

    def __str__(self):
        return "sort_reads()"

    @staticmethod
    def valid_path(readpath, groupfile):


        if not os.path.exists(readpath):
            raise FileNotFoundError("Read path is not found!")

        if not os.path.exists(groupfile):
            raise FileNotFoundError("Grouping file is not found!")




    @staticmethod
    def read_labels(readslist):
        read_pattern = re.compile("^lane1-.*-[ATCG]+?-[ATCG]+?-([a-zA-Z0-9-]+?)_")
        
        read_label_mapping = defaultdict(list)
        
        for i, v in enumerate(readslist):
                for m in [read_pattern.search(v)]:
                    if m:
                        read_label_mapping[re.sub("-|_", "", m.group(1))].append(v) 

        return read_label_mapping

    @staticmethod
    def matching_reads(label_to_match, mapping_labels):

        """
        mapping_labels - read pattern: reads
        returns a list of matched reads to that group
        """
        return mapping_labels[label_to_match]


    def grouping_reads(self, group_factor):

        # validate path
        try:
            Sort_reads.valid_path(self.readpath, self.groupfile)
        except FileNotFoundError as fnot:
            print(fnot)
            return "Input path are not found"


        meta = pd.read_csv(self.groupfile, sep="\t", names=["SampleID", "Site", "Location", "Flower", "Batch"], skiprows=1)

        # dict
        group_pattern = {}
        # grouping
        groupie = meta.groupby(group_factor)

        for i in groupie:
            group_pattern[i[0]] = list(i[1].SampleID.apply(lambda x: re.sub("-|_", "", x)))

        # sort reads

        reads_list = os.listdir(self.readpath)
        read_label_mapping = Sort_reads.read_labels(reads_list)

        for k, v in group_pattern.items():
            # v is the sample label in a list
            for j in v:

                read_matched = Sort_reads.matching_reads(j,read_label_mapping)

                # upacking read_matched (list of 2) read1 maybe reverse reads, here read1 and read2 do not have order
                read1 = os.path.join(self.readpath, read_matched[0])
                read2 = os.path.join(self.readpath, read_matched[1])

                self.grouped_reads[k].append(read1)
                self.grouped_reads[k].append(read2)

        return self


    def copy_reads(self, wkdir, subdirname):
        
        copy_dest = os.path.abspath(wkdir)

        #   reads
        #       |_ location1
        #       |_ location2
        #       |_ location3

        cur_dir = os.getcwd()
        os.chdir(copy_dest)

        if not os.path.exists(subdirname):
            os.makedirs(subdirname)

        base_path = os.path.abspath(subdirname)
        os.chdir(base_path)


        new_path = {}

        total_reads = len(os.listdir(self.readpath))
        print(f"Original total number of reads: {total_reads}\n")
        print(f"There are {len([ j for i in self.grouped_reads.values() for j in i])} reads in grouped files")


        

        for k, v in self.grouped_reads.items():
            # k is the group factor
            # v is the absolute path where the reads are stored in a list
            
            if not os.path.exists(k):
                os.makedirs(k)

            k_path = os.path.abspath(k)
            new_path[k] = k_path

            for i in tqdm(v, desc=f"{k} # reads - {len(v)} - progress"):

                shutil.copy2(i, k_path)

            


        os.chdir(cur_dir)

        print("Completed")

        return new_path

        
        
