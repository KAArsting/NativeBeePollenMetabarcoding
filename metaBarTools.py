import pandas as pd
import os
import re
import subprocess
import numpy as np
import shutil

"""[summary]
modified 7/25/2019
modified content: add matchedby parameter to manifest function so user can decide what match pattern to use
"""


def getPairedReads(mapping_index, readslist):
    """
    mapping_index: a dictionary-- folder for grouping: read label
    """
    all_reads_label = []
    for k, v in mapping_index.items():
        for i in v:
            all_reads_label.append(i)

    # find the matching pair in tuples
    def getKey(x):
        pattern = re.compile('(R[1-2]+)')
        y = pattern.search(x).group(0)
        return y


    reads_in_pair = []
    for match in all_reads_label:

        match_pattern = ".*("+match+").*"
        regex = re.compile(match_pattern)
        read_pair = [m.group() for read in readslist for m in [regex.search(read)] if m]

        if len(read_pair) >1:
            read_pair_sort = sorted(read_pair, key= lambda x: getKey(x))
            reads_in_pair.append(tuple(read_pair_sort))
        else:
            reads_in_pair.append(read_pair[0])

    return reads_in_pair



class metaBar_PreX:

    def __init__(self):
        self

    def metaBar_Copy(self, from_dir, to_dir, all_items = True):
        """
        copy from one to the other directories

        absolute path

        return to_dir
        """
        print("""
        ####################################################
        #########     metaBar_Copy                  ########
        ####################################################
        """)

        if not os.path.exists(from_dir):
            print("Copy From Path does not exist. Please check\n")
            return

        
        if not os.path.exists(to_dir):
            os.makedirs(to_dir)
            
        to_dir_path = os.path.abspath(to_dir)

        if all_items:
            print("\tCoping all items from {}\n".format(from_dir))
            allitem = from_dir + "/*"
            cmd = " ".join(["cp", allitem, to_dir_path])
            subprocess.run(cmd, shell=True, cwd=to_dir_path)

        return to_dir


    def metaBar_makeSubDir(self, uplevelname, *dir_sublist):
        """
        This function always makes sub directories based on the names in the list

        return path
        """
        print("""
        ####################################################
        #########          metaBar_makeSubDir     ##########
        ####################################################
        """)

        if not os.path.exists(uplevelname):
            try:
                os.makedirs(uplevelname)
            except os.error as e:
                print(e)
        
        uplevel_path = os.path.abspath(uplevelname)

        if len(dir_sublist) > 0:

            sublevels_path = list(map(lambda x: '/'.join([uplevel_path, x]), dir_sublist[0]))

            for i in sublevels_path:
                if not os.path.exists(i):
                    os.makedirs(i)
        else:
            sublevels_path = uplevel_path

        return sublevels_path


    def metaBar_LoadPlateSetup(self, platefile, primer_file, sheetname=0, colranges = [0,9], colnames = ['Plate', 'ID_No', 'Sample_ID', 'PCR_Conc', 'nmol_per_sample', 'Amount_of_Sample', 'Amount_of_Water', 'Well_No', 'Primer_set'], groupby =['Primer_set']):
        
        """
        This function take in plate setup csv files and output a mapping file to group reads, for example, by the primers used.

        platefile: a excel file or other kind
        primer_file: the file contains primers, format as below
        primer_name (corresponding to platefile):seq

        sheetname: sheetname to parse, integer
        colranges: select which columns you want
        colnames: parse the column names
        """
        print("""
        ####################################################
        #########     metaBar_LoadPlateSetup      ##########
        ####################################################
        """)

        # aggregate platefile
        useCol = list(range(colranges[0], colranges[-1]))
        plate = pd.read_excel(platefile, sheet_name=sheetname, usecols=useCol, names=colnames, dtype='object')

        # remove blank lines
        plate.dropna(axis=0, inplace=True)
        plate.reset_index(level=0, drop=True, inplace=True)
        # plate order associated with which reads they represent
        orig_index = list(plate.index.astype(dtype='int64'))
        # create a label that can match to read label
        pos_index = list(map(lambda x: "s00" + str(x+1) if x+1 <=9 else ("s0" + str(x+1) if x+1 <= 99 else "s" + str(x+1)), orig_index))

        plate["Read_label"] = pos_index

        # join the " " in Sample ID with "_"
        plate['Sample_ID'] = plate['Sample_ID'].apply(lambda x: "_".join(x.split()))
        # dealing with Primer_set &
        # 'and', ' and ', ' & ', '&'
        plate['Primer_set'] = plate['Primer_set'].apply(lambda x: "@".join(re.compile("_&_|_and_| & | and |&|and").split(x.replace(" ", ""))))


        # load primer file:
        primers = pd.read_table(primer_file, header = None, names=['oligos'])
        # get individual primers
        oligo,seqs = primers['oligos'].apply(lambda x:x.split(":")[0]),primers['oligos'].apply(lambda x:x.split(":")[1])
        # create primer and seq dictionary
        oligo_dic = {}
        for i in range(len(oligo)):
            oligo_dic[oligo[i]] = seqs[i]


        # output new plate setup
        subplate = plate.groupby(groupby)
        # the name should show reads labels (first-last) and groupby
        outfolder_names = []
        primer_F = []
        primer_R = []
        F_seq = []
        R_seq = []
        mapping_index = {}

        for key, _ in subplate:
            sub_df = subplate.get_group(key)
            # output folder name
            # label match
            label_match = []
            for i in sub_df['Read_label']:
                label_match.append(i)

            folder_name = '-'.join([sub_df['Read_label'].values[0],  sub_df['Read_label'].values[-1]])
            # key can be tuple if more than two
            if type(key) == str:
                key = (str(key),)

            key_list = [str(i) for i in key]
            key_join = '-'.join(key_list)
            folders = '_'.join([folder_name, key_join])
            outfolder_names.append(folders)
            mapping_index[folders] = label_match

            # get primer
            f_primer = list(map(lambda x: x.split('@')[0], sub_df['Primer_set'].unique()))
            r_primer = list(map(lambda x: x.split('@')[-1], sub_df['Primer_set'].unique()))

            primer_F.append(f_primer[0])
            primer_R.append(r_primer[0])
            F_seq.append(oligo_dic[f_primer[0]])
            R_seq.append(oligo_dic[r_primer[0]])

        map_out = pd.DataFrame({"Sample_Set":outfolder_names, "Forward_primer":primer_F, "Forward_seq":F_seq, "Reverse_primer":primer_R, "Reverse_seq":R_seq}, columns=['Sample_Set', 'Forward_primer', 'Forward_seq', 'Reverse_primer', 'Reverse_seq'])

        map_out.to_csv("NewPlateSetup.txt", sep="\t", header=True, index=False)


        return mapping_index

    def metaBar_Qiime2_Manifest(self, readpath, platefile, sheetname=0, colranges = [0,9], colnames = ['Plate', 'ID_No', 'Sample_ID', 'PCR_Conc', 'nmol_per_sample', 'Amount_of_Sample', 'Amount_of_Water', 'Well_No', 'Primer_set'], groupby =['Primer_set'], paired = True, matchby="index"):
        print("""
        ####################################################
        #########     metaBar_Qiime2_Manifest     ##########
        ####################################################
        """)
        # aggregate platefile
        useCol = list(range(colranges[0], colranges[-1]))
        plate = pd.read_excel(platefile, sheet_name=sheetname, usecols=useCol, names=colnames, dtype='object')

        # remove blank lines
        plate.dropna(axis=0, inplace=True)
        plate.reset_index(level=0, drop=True, inplace=True)
        # plate order associated with which reads they represent
        orig_index = list(plate.index.astype(dtype='int64'))
        # create a label that can match to read label
        pos_index = list(map(lambda x: "s00" + str(x+1) if x+1 <=9 else ("s0" + str(x+1) if x+1 <= 99 else "s" + str(x+1)), orig_index))

        plate["Read_label"] = pos_index

        # join the " " in Sample ID with "_"
        plate['Sample_ID'] = plate['Sample_ID'].apply(lambda x: "_".join(re.split(" |/", x)))
        # dealing with Primer_set &
        # 'and', ' and ', ' & ', '&'
        plate['Primer_set'] = plate['Primer_set'].apply(lambda x: "@".join(re.compile("_&_|_and_| & | and |&|and").split(x.replace(" ", ""))))


        subplate = plate.groupby(groupby)


        readsList = os.listdir(readpath)
        readspattern = [ re.sub("-|_", "", m.group(1)) for i in readsList for m in [re.search("^lane1-.*-[ATCG]+?-[ATCG]+?-([a-zA-Z0-9-]+?)_", i)] if m]



        if paired:

            for key, _ in subplate:

                group_df = subplate.get_group(key)

                if matchby == "index":
                    match_index = np.asarray(group_df['Read_label'])
                elif matchby == "sample":
                    match_index = list(map(lambda x: re.sub("-|_", "", x), np.asarray(group_df['Sample_ID'])))


                sample_names = np.asarray(group_df['Sample_ID'])

                sample_id = []
                absolut_filepath = []
                direction = []


                for i in range(len(match_index)):

                    if matchby == "index":
                        match_pattern = ".*("+match_index[i]+").*"
                        regex = re.compile(match_pattern)
                        read_pair = [m.group(0) for ind, read in enumerate(readsList) for m in [regex.search(read)] if m]

                    elif matchby == "sample":
                        match_pattern = "^"+match_index[i] + "$"
                        regex = re.compile(match_pattern)
                        read_pair = [ read for ind, read in enumerate(readsList) for m in [regex.search(readspattern[ind])] if m]

                    for j in read_pair:

                        sample_id.append(sample_names[i])
                        absolut_filepath.append(readpath+'/'+j)
                        if 'R1' in j:
                            direction.append('forward')
                        elif 'R2' in j:
                            direction.append('reverse')

                p_df = pd.DataFrame({'sample-id':sample_id, 'absolute-filepath':absolut_filepath, 'direction':direction}, columns=['sample-id', 'absolute-filepath', 'direction'])

                print("# of reads: {0}\n# of reads in the manifest file: {1}\n".format(len(readsList),len(p_df)))
                p_df.to_csv('{}_manifest.csv'.format(key), index=False, header = True)
        else:

            match_R1 = re.compile(".*_R1_.*")
            readsList_R1 = [ i for i in readsList if match_R1.match(i)]
            readspattern_R1 = [ re.sub("-|_", "", m.group(1)) for i in readsList_R1 for m in [re.search("^lane1-.*-[ATCG]+?-[ATCG]+?-([a-zA-Z0-9-]+?)_", i)] if m]

            for key, _ in subplate:

                group_df = subplate.get_group(key)

                if matchby == "index":
                    match_index = np.asarray(group_df['Read_label'])
                elif matchby == "sample":
                    match_index = list(map(lambda x: re.sub("-|_", "", x), np.asarray(group_df['Sample_ID'])))

                
                
                sample_names = np.asarray(group_df['Sample_ID'])

                sample_id = []
                absolut_filepath = []
                direction = []


                for i in range(len(match_index)):

                    if matchby == "index":
                        match_pattern = ".*("+match_index[i]+").*"
                        regex = re.compile(match_pattern)
                        read_pair = [m.group(0) for ind, read in enumerate(readsList_R1) for m in [regex.search(read)] if m]

                    elif matchby == "sample":
                        match_pattern = "^"+match_index[i] + "$"
                        regex = re.compile(match_pattern)
                        read_pair = [ read for ind, read in enumerate(readsList_R1) for m in [regex.search(readspattern_R1[ind])] if m]

                    for j in read_pair:

                        sample_id.append(sample_names[i])
                        absolut_filepath.append(readpath+'/'+j)
                        if 'R1' in j:
                            direction.append('forward')
                        elif 'R2' in j:
                            direction.append('reverse')

                p_df = pd.DataFrame({'sample-id':sample_id, 'absolute-filepath':absolut_filepath, 'direction':direction}, columns=['sample-id', 'absolute-filepath', 'direction'])
                
                print("# of reads: {0}\n# of reads in the manifest file: {1}\n".format(len(readsList_R1),len(p_df)))
                p_df.to_csv('{}_manifest.csv'.format(key), index=False, header = True)


        print("Manifest Completed")
        return "Completed"

                    

    def metaBar_Grouping(self, readpath, sublevels_path, mapping_index):
        """
        Use this to group reads based on their sequence ID
        read: path to the reads
        mapping_index: a dictionary contains folder name as the key, and match regex as values
        """

        print("""
        ####################################################
        #########          metaBar_Grouping      ###########
        ####################################################
        """)
        try:
            os.path.exists(readpath)
        except:
            print("Read path does not exists")
            return
        
        reads = os.listdir(readpath)

        if len(sublevels_path) > 0:
            paths = [(i.split('/')[-1], i) for i in sublevels_path]
        else:
            print("Folders to be copied to is empty")
            return


        for i in range(len(paths)):

            folder = paths[i][0]
            path = paths[i][-1]
            # get regex
            match_label = mapping_index[folder]

            for match in match_label:
                match_pattern = ".*("+match+").*"
                regex = re.compile(match_pattern)
                matched_reads = [m.group() for read in reads for m in [regex.search(read)] if m]

                for j in matched_reads:
                    print("Coping reads to your folder: {}".format(folder))
                    cp_path = '/'.join([readpath, j])
                    cmd = ' '.join(["cp", cp_path, path])
                    subprocess.run(cmd, shell=True)

        print("\nGrouping Completed!\n")
        return


    def metaBar_TrimReads(self,mapping_index, readpath, Trim2Folder, TrimSoftware="$TRIMMO", ReadMode="PE", threads = 1, **kwargs):
        """
        To cp all reads from the path and trimmed them use trimmomatic
        readpath: path to reads
        Trim2Folder: Folder to contain the trimmed reads (paired and unpaired as sub directoryies)
        TrimSoftware: absolute path of trimmo, you can put alias in your path, such as $TRIMMO

        **kwargs: use for trimmomatic
        """
        print("""
        ####################################################
        #########          Trimmomatic      ################
        ####################################################
        """)

        try:
            os.path.exists(readpath)
        except:
            print("Raw reads path not existed")
            return

        try:
            os.path.exists(Trim2Folder[0])
        except:
            print("Folder for trim reads not existed")
            return




        if len(kwargs) > 0:


            reads = os.listdir(readpath)
            # with trange(len(reads)) as t:
            #     for i in t:
            #         cp_link = readpath + "/" + reads[i]
            #         cmd = " ".join(["cp", cp_link, Trim2Folder])
            #         subprocess.run(cmd, shell=True)


            print("Starting Triming reads ...... \n")


            # get paired reads in tuple if PE, list if SE
            pairs = getPairedReads(mapping_index, reads)
            
            if ReadMode == "PE":

                trim_params = []
                for k, v in kwargs.items():
                    param = ':'.join([k, str(v)])
                    trim_params.append(param)

                for i in pairs:
                    read_f = readpath + '/' + i[0]
                    read_r = readpath + '/' + i[1]

                    paired_folder = Trim2Folder[0]
                    unpaired_folder = Trim2Folder[1]
                    fout_paired = paired_folder+"/"+"trimmed_PE_"+i[0]
                    fout_unpaired = unpaired_folder+"/"+"trimmed_unPE_"+i[0]
                    rout_paired = paired_folder+"/"+"trimmed_PE_"+i[1]
                    rout_unpaired = unpaired_folder + "/"+"trimed_unPE_"+i[1]

                    trim_cmd = ["java -jar", TrimSoftware, ReadMode, "-threads", str(threads), read_f, read_r, fout_paired, fout_unpaired, rout_paired, rout_unpaired]
                    trim_cmd.extend(trim_params)
                    trim_cmd_all = " ".join(trim_cmd)
                    print(trim_cmd_all)

                    try:
                        subprocess.run(trim_cmd_all, shell=True, cwd=Trim2Folder)
                    except subprocess.SubprocessError as e:
                        print(e)

                return paired_folder
                
            elif ReadMode == "SE":

                for i in pairs:
                    read_se = readpath + "/" + i
                    out = Trim2Folder+"/"+"trimmed_SE_" + i
                    trim_cmd = ["java -jar", TrimSoftware, ReadMode, "-threads", str(threads), read_se, out]
                    trim_cmd.extend(trim_params)
                    trim_cmd_all = " ".join(trim_cmd)
                    try:
                        subprocess.run(trim_cmd_all, shell=True, cwd=Trim2Folder)
                    except subprocess.SubprocessError as e:
                        print(e)

                return out

            print("Trimming Reads Completed!\n")

        else:
            print("Please Enter Parameters For Trimmomatic")
            return


class metaBar_Vsearch:
    
    def __init__(self):

        self.flag = []
        #self.readpath = readpath
        print("Checking if vsearch is installed in your path\n")

        if shutil.which('vsearch') == None:
            print("Vsearch is not installed in the path")
            return
        else:
            print("Vsearch is installed, proceed to the use the pipeline >>>")


    def __str__(self):
        steps = [ "Step{}: {}".format(i+1, self.flag[i]) for i in range(len(self.flag))]
        pipeline = "--->".join(steps)
        return pipeline


    def vsearch_derep(self, input_path, wkpath, input_param='--derep_fulllength', output_param='--output', **kwargs):
        """
        dereplicate the reads in the path

        input_path: path to the input reads
        wkpath: the folder to hold the results

        input must be fasta file
        """

        print("""
        ####################################################
        #########      Vsearch Dereplications      #########
        ####################################################
        """)
        print("\n")
        print("""
        vsearch (--derep_fulllength | --derep_prefix) fastafile (--output | --uc) outputfile [options]
        
        you don't have to specify relabel string, just put down --relabel
        """)

        readsList = os.listdir(input_path)

        if len(readsList) == 0:
            print("Empty input file. Exit Program")
            return 

        if not readsList[0].endswith('fasta'):
            print("Input reads should be fasta format, please check your pipeline\n")
            return
        else:
            for i in readsList:

                if len(kwargs) <= 0:
                    output_name = "derep_" + "".join(i.split(".")[0])
                    read = input_path + "/" + i

                    if output_param == "--output":
                        outfile = output_name+".fasta"
                    elif output_param == "--uc":
                        outfile = output_name + ".uc"

                    cmd = ["vsearch", input_param, read, output_param, outfile]
                    input_cmds = " ".join(cmd)

                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

                else:

                    output_name = "derep_" + "".join(i.split(".")[0])
                    read = input_path + "/" + i

                    if output_param == "--output":
                        outfile = output_name+".fasta"
                    elif output_param == "--uc":
                        outfile = output_name + ".uc"

                    params = []
                    for k, v in kwargs.items():

                        cmd = "--" + k
                        par = str(v)

                        if cmd == input_param or cmd == output_param:
                            print("Repeated Parameters. Please Check input parameters")
                            return
                        
                        if par == 'None':
                            if cmd == "--uc":
                                outfile2 = output_name + ".uc" 
                                params.append(" ".join([cmd, outfile2]))

                            elif cmd == "--output":
                                outfile2 = output_name + ".fasta" 
                                params.append(" ".join([cmd, outfile2]))
                            
                            elif cmd == "--relabel":
                                # get relabel string based on sequence Identifier s001, s002, ...
                                pattern = "(s[0-9]+)"
                                regex = re.compile(pattern)
                                m = [l.group(0) for l in [regex.search(i)] if l ]
                                string = "".join(m)

                                params.append(" ".join([cmd, "OTU_"+string]))

                            else:
                                params.append(" ".join([cmd]))
                        else:
                            params.append(" ".join([cmd, par]))


                    cmds = ["vsearch", input_param, read, output_param, outfile]
                    cmds.extend(params)
                    input_cmds = " ".join(cmds)
                    print(input_cmds)
                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

        print("Completed Dereplication\n")
        self.flag.append("Dereplication")
        return wkpath

    def vsearch_merge(self, mapping_index, input_path, wkpath, output_param = '--fastaout', **kwargs):

        """
        merge paire end reads
        """
        print("""
        ####################################################
        #########      Vsearch Mergepaires        ##########
        ####################################################
        """)

        print("""vsearch --fastq_mergepairs fastqfile --reverse fastqfile (--fastaout | --fastqout | --fastaout_not-
merged_fwd | --fastaout_notmerged_rev | --fastqout_notmerged_fwd | --fastqout_notmerged_rev |
--eetabbedout) outputfile [options]
""")
        print("For output, you don't need to specify the name, just do fastaout='None' or fastaout_not-merged_fwd='None'")

        # get paired reads
        readslist = os.listdir(input_path)
        if len(readslist) == 0:
            print("Empty input file. Exit Program")
            return

        paired_reads = getPairedReads(mapping_index, readslist)
 
        # separated reads by R1 and R2
        R1_reads = [paired_reads[i][0] for i in range(len(paired_reads))]
        R2_reads = [paired_reads[i][1] for i in range(len(paired_reads))]


        for i in range(len(R1_reads)):
            r1 = R1_reads[i]
            r2 = R2_reads[i]

            output_name = "merged_" + "".join(r1.split(".")[0])
            print("\n")
            print("Forward Read: {}".format(r1))
            print("Reverse Read: {}".format(r2))

            if len(kwargs) <= 0:
                if output_param == "--fastaout":
                    outfile = output_name + ".fasta"
                elif output_param == "--fastqout":
                    outfile = output_name + ".fastq"         

                cmds = ["vsearch --fastq_mergepairs", input_path+"/"+r1, "--reverse", input_path+"/"+r2, output_param, outfile]
                input_cmds = " ".join(cmds)
                try:
                    subprocess.run(input_cmds, shell=True, cwd=wkpath)
                except subprocess.SubprocessError as e:
                    print(e)

            elif len(kwargs) > 0:

                if output_param == "--fastaout":
                    outfile = output_name + ".fasta"
                elif output_param == "--fastqout":
                    outfile = output_name + ".fastq" 

                params = []
                for k, v in kwargs.items():

                    cmd = "--" + k
                    par = str(v)

                    if cmd == output_param or cmd == '--fastq_mergepairs':
                        print("Repeated output parameters. Please Check")
                        return

                    if par == 'None':
                        if cmd == "--fastaout":
                            outfile2 = output_name + ".fasta"
                        elif cmd == "--fastqout":
                            outfile2 = output_name + ".fastq"
                        elif cmd == "--fastaout_not-merged_fwd":
                            outifle2 = "Not-fwd"+output_name + ".fasta"
                        elif cmd == "--fastaout_not-merged_rev":
                            outifle2 = "Not-rev"+output_name + ".fasta"
                        elif cmd == "--fastqout_not-merged_fwd":
                            outifle2 = "Not-fwd"+output_name + ".fastq"
                        elif cmd == "--fastqout_not-merged_rev":
                            outifle2 = "Not-rev"+output_name + ".fastq"
                        else:
                            outfile2 = ""

                        if outfile2 == "":
                            params.append(" ".join([cmd]))
                        else:
                            params.append(" ".join([cmd, outfile2]))
                    else:
                        params.append(" ".join([cmd, par]))
    
                cmds = ["vsearch --fastq_mergepairs", input_path+"/"+r1, "--reverse", input_path+"/"+r2, output_param, outfile]
                cmds.extend(params)
                input_cmds = " ".join(cmds)
                try:
                    subprocess.run(input_cmds, shell=True, cwd=wkpath)
                except subprocess.SubprocessError as e:
                    print(e)


        print("Merging Completed!")
        self.flag.append("Merge")
        return wkpath

    def vsearch_quality(self, input_path, wkpath, input_param = "--fastq_eestats", **kwargs):

        """
        input fastq file
        """

        print("""
        ####################################################
        #########   Vsearch eestats check quality ##########
        ####################################################
        """)

        readsList = os.listdir(input_path)
        if len(readsList) == 0:
            print("Empty Input File. Exit Program")
            return

        if not readsList[0].endswith(('fastq', 'fastq.gz')):
            print("Input reads should be fasta format, please check your pipeline\n")
            return

        else:
            for i in readsList:
                print("Calculate quality statistics\n")
                if len(kwargs) <= 0:
                    output_name = "stats_" + "".join(i.split(".")[0]) + ".stats"

                    read = input_path + '/' + i
                    cmd = ["vsearch", input_param, read, "--output", output_name]
                    input_cmds = " ".join(cmd)

                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

                else:
                    params = []
                    for k, v in kwargs.items():

                        cmd = "--" + k
                        par = str(v)
                        if par == 'None':
                            params.append(cmd)
                        else:
                            params.append(" ".join([cmd, par]))

                    output_name = "stats_" + "".join(i.split(".")[0]) + ".stats"
                    read = input_path + '/' + i
                    cmd = ["vsearch", input_param, read, "--output", output_name]
                    input_cmds = " ".join(cmd)

                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

        print("Completed checking quality statistics\n")
        self.flag.append("quality_stats")
        return wkpath


    def vsearch_filter(self, input_path, wkpath, output_param = '--fastaout', **kwargs):
        """
        input_path: input fastq file
        output_param: defines why type of output you want to
        """
        print("""
        ####################################################
        #########   Vsearch filter fastq          ##########
        ####################################################
        """)

        print("""vsearch --fastq_filter fastqfile (--fastaout | --fastaout_discarded | --fastqout | --fastqout_discarded)
outputfile [options]\n""")

        readsList = os.listdir(input_path)
        if len(readsList) == 0:
            print("Empty Input File. Exit Program")
            return

        if not readsList[0].endswith(('fastq', 'fastq.gz')):
            print("Input reads should be fastq format, please check your pipeline\n")
            return


        else:
            for i in readsList:


                output_name = "filtered_" + "".join(i.split(".")[0])
                if 'fasta' in output_param:
                    outfile = output_name+".fasta"
                elif 'fastq' in output_param:
                    outfile = output_name+".fastq"

                read = input_path + '/' + i


                if len(kwargs) <= 0:
                    cmd = ["vsearch", "--fastq_filter", read, output_param, outfile]
                    input_cmds = " ".join(cmd)

                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

                else:
                    params = []
                    for k, v in kwargs.items():

                        cmd = "--" + k
                        par = str(v)

                        if cmd == output_param:
                            print("Repeated Parameters. Please Check.")
                            return

                        if par == 'None':
                            if cmd == "--fastaout":
                                outfile2 = output_name + ".fasta"
                            elif cmd == "--fastaout_discarded":
                                outfile2 = "Discarded_"+output_name+".fasta"
                            elif cmd == "--fastqout":
                                outfile2 = output_name+".fastq"
                            elif cmd == "--fastqout_discarded":
                                outfile2 = "Discarded_"+output_name+".fastq"
                            else:
                                outfile2 = ""
                            
                            if outfile2 == "":
                                params.append(cmd)
                            else:
                                params.append(" ".join([cmd, outfile2]))
                        else:
                            params.append(" ".join([cmd, par]))

                    if 'fasta' in output_param:
                        outfile = output_name+".fasta"
                    elif 'fastq' in output_param:
                        outfile = output_name+".fastq"

                    cmd = ["vsearch", "--fastq_filter", read, output_param, outfile]
                    cmd.extend(params)
                    input_cmds = " ".join(cmd)

                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

        print("Completed Filtering\n")
        self.flag.append("Filtering")
        return wkpath

    def vsearch_chimera_uchime(self, input_path, wkpath, uchime = '--uchime3_denovo', uchime_out = "--nonchimeras", **kwargs):
    
        print("""
        ####################################################
        ######### Vsearch chimera  uchime3_denovo ##########
        ####################################################
        """)

        print("""
        vsearch (--uchime_denovo | --uchime2_denovo | --uchime3_denovo) fastafile (--chimeras |
        --nonchimeras | --uchimealns | --uchimeout) outputfile [options]\n
        vsearch --uchime_ref fastafile (--chimeras | --nonchimeras | --uchimealns | --uchimeout) outputfile
        --db fastafile [options]
        """)

        print("You don't have to specify fileanme for output parameters, simply do --chimera\n")

        readsList = os.listdir(input_path)
        if len(readsList) == 0:
            print("Empty Input File. Exit Program")
            return


        if not readsList[0].endswith('fasta'):
            print("Input reads should be fasta format, please check your pipeline\n")
            return

        else:
            for i in readsList:

                output_name = "".join(i.split(".")[0])
                read = input_path + '/' + i

                if len(kwargs) <= 0:
                    
                    if uchime_out == "--nonchimeras":
                        outfile = "NonChimeras_"+output_name+".fasta"
                    elif uchime_out == "--chimeras":
                        outfile = "Chimeras_"+output_name+".fasta"
                    elif uchime_out == "--uchimealns":
                        outfile = "Chimera_"+output_name+".align"
                    elif uchime_out == "--uchimeout":
                        outfile = "Uchimeout_"+output_name+".tab"

                    cmd = ["vsearch", uchime, read, uchime_out, outfile]
                    input_cmds = " ".join(cmd)

                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

                else:
                    params = []
                    for k, v in kwargs.items():

                        cmd = "--" + k
                        par = str(v)

                        if cmd == uchime_out or cmd == uchime:
                            print("Repeated Parameters. Please Check")
                            return


                        if par == "None":
                            if cmd == "--nonchimeras":
                                outfile2 = "NonChimeras_"+output_name+".fasta"
                            elif cmd == "--chimeras":
                                outfile2 = "Chimeras_"+output_name+".fasta"
                            elif cmd == "--uchimealns":
                                outfile2 = "Chimera_"+output_name+".align"
                            elif cmd == "--uchimeout":
                                outfile2 = "Uchimeout_"+output_name+".tab"
                            else:
                                outfile2 = ""

                            if outfile2 == "":
                                params.append(cmd)
                            else:
                                params.append(" ".join([cmd, outfile2]))
                        else:
                            params.append(" ".join([cmd, par]))

                    if uchime_out == "--nonchimeras":
                        outfile = "NonChimeras_"+output_name+".fasta"
                    elif uchime_out == "--chimeras":
                        outfile = "Chimeras_"+output_name+".fasta"
                    elif uchime_out == "--uchimealns":
                        outfile = "Chimera_"+output_name+".align"
                    elif uchime_out == "--uchimeout":
                        outfile = "Uchimeout_"+output_name+".tab"

                    cmds = ["vsearch", uchime, read, uchime_out, outfile]
                    cmds.extend(params)
                    input_cmds = " ".join(cmds)

                    try:
                        subprocess.run(input_cmds, shell=True, cwd = wkpath)
                    except subprocess.SubprocessError as e:
                        print(e)

        print("Completed Chimera Detection\n")
        return wkpath



    def vsearch_cluster(self, input_path, wkpath, cluster_method = '--cluster_fast', output_param = '--otutabout', id = 0.95, **kwargs):
        """
        cluster fasta

        Clustering:
        vsearch (--cluster_fast | --cluster_size | --cluster_smallmem | --cluster_unoise) fastafile (--alnout |
        --biomout | --blast6out | --centroids | --clusters | --mothur_shared_out | --msaout | --otutabout |
        --profile | --samout | --uc | --userout) outputfile --id real [options]

        """


        print("""
        ####################################################
        #########       Vsearch cluster           ##########
        ####################################################
        """)

        print("\n")
        print("""
        vsearch (--cluster_fast | --cluster_size | --cluster_smallmem | --cluster_unoise) fastafile (--alnout |
        --biomout | --blast6out | --centroids | --clusters | --mothur_shared_out | --msaout | --otutabout |
        --profile | --samout | --uc | --userout) outputfile --id real [options]
        """)

        fasta = os.listdir(input_path)

        if len(fasta) == 0:
            print("Empty Input File. Exit Program")
            return

        if not fasta[0].endswith('fasta'):
            print("Input reads should be fasta format, please check your pipeline\n")
            return

        for i in fasta:

            output_name = "clustered_" + "".join(i.split(".")[0])
            read = input_path + '/' + i


            if len(kwargs) <= 0:

                if output_param == "--centroids":
                    outfile = "Centroids_"+output_name+".fasta"
                elif output_param == "--consout":
                    outfile = "Consensus_"+output_name+".fasta"
                elif output_param == "--biomout":
                    outfile = output_name+".biom"
                elif output_param == "--otutabout":
                    outfile = output_name+".otutab"
                elif output_param == "--uc":
                    outfile = output_name+".uc"
                elif output_param == "--msaout":
                    outfile = output_name+".msaout"
                elif output_param == "--mothur_shared_out":
                    outfile = output_name+".mothur_shared_out"
                elif output_param == "--blast6out":
                    outfile = output_name+".blast6out"
                elif output_param == "--clusters":
                    # get the string pattern. Use read identifier
                    pattern = "(s[0-9]+)"
                    regex = re.compile(pattern)
                    m = [l.group(0) for l in [regex.search(i)] if l ]
                    string = "".join(m)

                    outfile = string


                cmd = ["vsearch", cluster_method, read, output_param, outfile, "--id", str(id)]
                input_cmds = " ".join(cmd)


                try:
                    subprocess.run(input_cmds, shell=True, cwd = wkpath)
                except subprocess.SubprocessError as e:
                    print(e)

            else:
                params = []
                for k, v in kwargs.items():

                    cmd = "--" + k
                    par = str(v)

                    if cmd == cluster_method or cmd == output_param:
                        print("Repeated Parameters. Please Check")
                        return

                    if par == 'None':
                        if cmd == "--centroids":
                            outfile2 = "Centroids_"+output_name+".fasta"
                        elif cmd == "--consout":
                            outfile2 = "Consensus_"+output_name+".fasta"
                        elif cmd == "--biomout":
                            outfile2 = output_name+".biom"
                        elif cmd == "--otutabout":
                            outfile2 = output_name+".otutab"
                        elif cmd == "--uc":
                            outfile2 = output_name+".uc"
                        elif cmd == "--msaout":
                            outfile2 = output_name+".msaout"
                        elif cmd == "--mothur_shared_out":
                            outfile2 = output_name+".mothur_shared_out"
                        elif cmd == "--blast6out":
                            outfile2 = output_name+".blast6out"
                        elif cmd == "--clusters":
                            # get the string pattern. Use read identifier
                            pattern = "(s[0-9]+)"
                            regex = re.compile(pattern)
                            m = [l.group(0) for l in [regex.search(i)] if l ]
                            string = "".join(m)

                            outfile2 = string
                        else:
                            outfile2 = ""
                        
                        if outfile2 == "":
                            params.append(cmd)
                        else:
                            params.append(" ".join([cmd, outfile2]))
                    else:
                        params.append(" ".join([cmd, par]))

                if output_param == "--centroids":
                    outfile = "Centroids_"+output_name+".fasta"
                elif output_param == "--consout":
                    outfile = "Consensus_"+output_name+".fasta"
                elif output_param == "--biomout":
                    outfile = output_name+".biom"
                elif output_param == "--otutabout":
                    outfile = output_name+".otutab"
                elif output_param == "--uc":
                    outfile = output_name+".uc"
                elif output_param == "--msaout":
                    outfile = output_name+".msaout"
                elif output_param == "--mothur_shared_out":
                    outfile = output_name+".mothur_shared_out"
                elif output_param == "--blast6out":
                    outfile = output_name+".blast6out"
                elif output_param == "--clusters":
                    # get the string pattern. Use read identifier
                    pattern = "(s[0-9]+)"
                    regex = re.compile(pattern)
                    m = [l.group(0) for l in [regex.search(i)] if l ]
                    string = "".join(m)
                    outfile = string


                cmds = ["vsearch", cluster_method, read, output_param, outfile, "--id", str(id)]
                cmds.extend(params)
                input_cmds = " ".join(cmds)

                try:
                    subprocess.run(input_cmds, shell=True, cwd = wkpath)
                except subprocess.SubprocessError as e:
                    print(e)

        print("Completed Clustering\n")
        self.flag.append("Clustering")
        return wkpath

    def vsearch_search(self, input_path, wkpath, db, input_param="--usearch_global", output_param="--otutabout", id = 0.97, **kwargs):
        """
        can be used to generate OTU table
        """

        print("""
        ####################################################
        #########       Vsearch usearch_global    ##########
        ####################################################
        """)
        print("\n")
        print("""
        vsearch --usearch_global fastafile --db fastafile (--alnout | --biomout | --blast6out |
        --mothur_shared_out | --otutabout | --samout | --uc | --userout) outputfile --id real [options]
        """)


        fasta = os.listdir(input_path)

        if len(fasta) == 0:
            print("Empty Input File. Exit Program")
            return

        if not fasta[0].endswith('fasta'):
            print("Input reads should be fasta format, please check your pipeline\n")
            return

        for i in fasta:
            output_name = "".join(i.split(".")[0])
            read = input_path + '/' + i

            if len(kwargs) <=0:
                if output_param == "--biomout":
                    outfile = output_name+".biom"
                elif output_param == "--blast6out":
                    outfile = output_name+".biom"
                elif output_param == "--otutabout":
                    outfile = output_name+".otutab"
                elif output_param == "--uc":
                    outfile = output_name+".uc"
                
                cmd = ["vsearch", input_param, read, '--db', db, output_param, outfile, "--id", str(id)]
                input_cmds = " ".join(cmd)


                try:
                    subprocess.run(input_cmds, shell=True, cwd = wkpath)
                except subprocess.SubprocessError as e:
                    print(e)

        
            else:
                params = []
                for k, v in kwargs.items():

                    cmd = "--" + k
                    par = str(v)

                    if cmd == input_param or cmd == output_param:
                        print("Repeated Parameters. Please Check")
                        return

                    if par == 'None':
                        if cmd == "--alnout":
                            outfile = output_name+".alnout"
                        elif cmd == "--biomout":
                            outfile = output_name+".biom"
                        elif cmd == "--blast6out":
                            outfile = output_name+".blast6out"
                        elif cmd == "--otutabout":
                            outfile = output_name+".otutab"
                        elif cmd == "--uc":
                            outfile = output_name+".uc"
                        elif cmd == "--samout":
                            outfile = output_name+".samout"
                        elif cmd == "--mothur_shared_out":
                            outfile = output_name+".mothur_shared_out"
                        else:
                            outfile = ""
                        
                        if outfile == "":
                            params.append(cmd)
                        else:
                            params.append(" ".join([cmd, outfile]))
                    else:
                        params.append(" ".join([cmd, par]))

                        
                if output_param == "--biomout":
                    outfile = output_name+".biom"
                elif output_param == "--blast6out":
                    outfile = output_name+".biom"
                elif output_param == "--otutabout":
                    outfile = output_name+".otutab"
                elif output_param == "--uc":
                    outfile = output_name+".uc"

                cmds = ["vsearch", input_param, read, "--db", db, output_param, outfile, "--id", str(id)]
                cmds.extend(params)
                input_cmds = " ".join(cmds)

                try:
                    subprocess.run(input_cmds, shell=True, cwd = wkpath)
                except subprocess.SubprocessError as e:
                    print(e)
        
        print("Usearch Global Completed")
        self.flag.append("Usearch_Global")
        return wkpath

    def vsearch_taxon_assign(self, input_path, wkpath, db, cutoff = 0.97):
        """
        assign taxonomy
        """
        print("""
        ####################################################
        #########  Vsearch taxonomy classification  ########
        ####################################################
        """)
        print("\n")
        print("""
        The reference database must contain taxonomic information in the header of each sequence in the
        form of a string starting with ";tax=" and followed by a comma-separated list of up to eight taxo-
        nomic identifiers. Each taxonomic identifier must start with an indication of the rank by one of the
        letters d (for domain) k (kingdom), p (phylum), c (class), o (order), f (family), g (genus), or s
        (species). The letter is followed by a colon (:) and the name of that rank. Commas and semicolons
        are not allowed in the name of the rank.
        Example:
        ">X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacte-
        ria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli".
        """)
        seq = os.listdir(input_path)
        if len(seq) == 0:
            print("Empty Input File. Exit Program")
            return

        for i in seq:
            output_name = "clustered_" + "".join(i.split(".")[0])
            read = input_path + '/' + i

            outfile = output_name+'.taxout'

            cmd = ['vsearch --sintax', read, '--db', db, '--sintax_cutoff', str(cutoff), '--tabbedout', outfile]
            input_cmds = " ".join(cmd)

            try:
                subprocess.run(input_cmds, shell=True, cwd=wkpath)
            except subprocess.SubprocessError as e:
                print(e)
                return

        print("Taxonomy Classification Completed\n")
        self.flag.append("Taxonomy_classification")
        return wkpath
