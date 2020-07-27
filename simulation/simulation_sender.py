#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
Generates simulation and sends XICRA command
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import pandas as pd
from collections import defaultdict
import argparse
import numpy as np
from termcolor import colored

## import my modules
from XICRA.scripts import functions
from XICRA.scripts import reads2tabular
from XICRA.modules import prep
from XICRA.modules import join
from XICRA.modules import miRNA

## paired - end y luego hacer R1 y R2
## HS25 HiSeq 2500
## cobertura: 10

#########################################################
def NGS_simulator(name, abs_folder, seqSys_list, type_reads, fcov_list, fasta, 
                  threads_given, debug, art_illumina_bin, seqtk_bin, database_folder, send_XICRA):

    for profile in seqSys_list:
        print (" + Simulating reads for profile: " + profile)
        ## check how many profiles to use
        if len(seqSys_list) == 1:
            profile_path = abs_folder
        else:
            ## create folder
            profile_path = functions.create_subfolder(profile, abs_folder)

        for reads in type_reads:
            print ("  + Reads type: " + reads)
            ## check how many profiles to use
            if len(type_reads) == 1:
                reads_path = profile_path
            else:
                reads_path = functions.create_subfolder(reads, profile_path)

            for fcov in fcov_list:
                print ("   + Coverage x" + fcov)
                ## check how many profiles to use
                if len(fcov_list) == 1:
                    coverage_path = reads_path
                else:
                    coverage_path = functions.create_subfolder("x" + fcov, reads_path)

                ## simulate
                print ("    ***** Simulate NGS reads *****")
                print ("    ** Reads: " + reads)
                print ("    ** Built-in profile: " + profile)
                print ("    ** Read length: Determined by isomiR length")
                print ("    ** Fold coverage: " + fcov)
                print ("\n")
                
                ##
                tmp_path = functions.create_subfolder('tmp_fasta', coverage_path)
                tmp_fastq = functions.create_subfolder('tmp_fastq', coverage_path)
                
                # split isomiRs by length into multiple fasta files
                fasta_dict = process_fasta_length(fasta, tmp_path, debug)

                ## simulate reads for each length
                for int_len, fasta_file_len in fasta_dict.items():
                    str_len = str(int_len)
                    print ("\n    ** Simulate reads:")
                    print ("    ** Sequence fasta file: " + fasta_file_len)
                    print ("    ** Length: " + str_len)
                    
                    ## add random seed? --rndSeed for exact simulations?
                    ## outfile
                    outfile_name = name + '_P-' + profile + '_T-' + reads + '_x-' + fcov + '_L-' + str(int_len) + '_R'
                    outfile_path = os.path.join(tmp_fastq, outfile_name)
                    
                    ## command
                    art_illumina_cmd = "%s --rndSeed 123456789 -sam -ss %s -i %s -l %s -c %s -o %s" %(art_illumina_bin, profile, fasta_file_len, str_len, fcov, outfile_path)
                    
                    if (reads == 'PE'):
                        art_illumina_cmd = art_illumina_cmd + " -p -m 50 -s 5"
                    
                    
                    code = functions.system_call(art_illumina_cmd)
                    if not code:
                        print ("** ERROR: Some error happened during ART simulation")
                        exit()
                    
                    ## filter R1
                    fastq_R2_ids = discard_revcomp(outfile_path, reads)
                    
                    if (reads == 'PE'):
                        ## adjust R2
                        R2_fastq = outfile_path + '2.fq'
                        R2_filter = outfile_path + 'filter_R2.fq'
                        with open(R2_filter, 'w') as outfh:
                             with open(R2_fastq, 'r') as fh:
                                 lines = []
                                 for line in fh:
                                     lines.append(line.rstrip())
                                     if len(lines) == 4:
                                          record = reads2tabular.process_fastq(lines)
                                          #sys.stderr.write("Record: %s\n" % (str(record)))
                                          lines = []                                 
                                          if record['name'] in fastq_R2_ids.keys():
                                              outfh.write("%s\n%s\n%s\n%s\n" % (record['name'], record['sequence'], record['optional'], record['quality']))
                                                
                ## merge reads all lengths
                reads_path = functions.create_subfolder('reads', coverage_path)
                if (reads == 'PE'):
                    R1_all_reads = functions.retrieve_matching_files(tmp_fastq, "filter_R1.fq")
                    R1_reads = os.path.join(reads_path, name + '_R1.fq')
                else:
                    R1_all_reads = functions.retrieve_matching_files(tmp_fastq, "filter.fq")
                    R1_reads = os.path.join(reads_path, name + '.fq')
                    
                functions.merge_files(R1_reads, R1_all_reads)
                
                if (reads == 'PE'):
                    ## concat all reads 
                    R2_all_reads_tmp = functions.retrieve_matching_files(tmp_fastq, "filter_R2.fq")
                    R2_reads = os.path.join(reads_path, name + '_R2.fq')
                    functions.merge_files(R2_reads, R2_all_reads_tmp)

                ## 
                if send_XICRA:
                    if (reads == 'PE'):
                        call_XICRA_PE(coverage_path, reads_path, name, threads_given, debug, database_folder, seqtk_bin, R2_reads)
                        ## print time stamp
                        filename_stamp = coverage_path + '/.success'
                        functions.print_time_stamp(filename_stamp)
                    else:
                        call_XICRA_SE(coverage_path, reads_path, name, threads_given, debug, database_folder, seqtk_bin, R2_reads)
                        ## print time stamp
                        filename_stamp = coverage_path + '/.success'
                        functions.print_time_stamp(filename_stamp)
                else:
                    print ("+ Simulation for is ready in folder: ")
                    print (coverage_path)
    
    return()

##############
def discard_revcomp(outfile_path, reads):
    ##### Remove non 5'-3' simulated reads
    ## use art illumina aln file generated for R1
    if (reads == 'PE'):
        aln_file_R1 = outfile_path + '1.aln'
    else:
        aln_file_R1 = outfile_path + '.aln'
                    
    
    ## read aln file
    freq_fasta = defaultdict(int)
    fastq_dict = defaultdict(int)
    with open(aln_file_R1, 'r') as fh:
        lines = []
        for line in fh:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            if line.startswith('>'):
                line_list = line.rstrip().split('\t')
                if line_list[3] == '+':
                    ID=line_list[1]
                    lines.append(ID[:-2])
                continue
            else:
                if len(lines) == 1:
                    lines.append(line.rstrip())
                                      
            if len(lines) == 2:
                record = reads2tabular.process_fasta(lines)
                ##sys.stderr.write("Record: %s\n" % (str(record)))
                lines = []
                    
                ## add sequences & count
                freq_fasta[record['sequence']] += 1
    
    if (reads == 'PE'):
        ## read R1 fastq file
        fastq_file = outfile_path + '1.fq'
        out_file = outfile_path + '_filter_R1.fq'
    else:
        fastq_file = outfile_path + '.fq'
        out_file = outfile_path + '_filter.fq'
        
        
    ## print in file
    with open(out_file,'w') as file:
        with open(fastq_file, 'r') as fh:
            lines = []
            for line in fh:
                lines.append(line.rstrip())
                if len(lines) == 4:
                    record = reads2tabular.process_fastq(lines)
                    #sys.stderr.write("Record: %s\n" % (str(record)))
                    lines = []
                    fastq_ID = record['name'].replace('/1', '/2')
                    
                    if record['sequence'] in freq_fasta.keys():
                        file.write("%s\n%s\n%s\n%s\n" % (record['name'], record['sequence'], record['optional'], record['quality']))
                        fastq_dict[fastq_ID] += 1
    
    file.close()
    return (fastq_dict)

###################
def call_XICRA_SE(folder_path, reads_path, name, threads_given, debug_bool, database_folder, seqtk_bin, R2_reads):
    
    ## debugging messages
    if debug_bool:
        print ("\n********* XICRA SE analysis *********\n")
    
    ## send XICRA command
    ## create argparse with arguments provided to call XICRA miRNA
    output_folder_XICRA = os.path.join(folder_path, 'analysis_R1')
    
    software_miRNA = ["sRNAbench", "optimir", "miraligner"]
    XICRA_options_miRNA = argparse.Namespace(input=reads_path, output_folder=output_folder_XICRA,
                                            single_end=True, batch=False, in_sample=False, noTrim=True, 
                                            ex_sample=False, detached=False, include_lane=False, 
                                            include_all=False, threads=threads_given,
                                            soft_name=software_miRNA, species='hsa',
                                            database=database_folder, miRNA_gff=False, hairpinFasta=False, 
                                            matureFasta=False, miRBase_str=False, 
                                            help_format=False,  help_project=False, help_miRNA=False, debug=debug_bool)
    miRNA.run_miRNA(XICRA_options_miRNA)
    return()

###################
def call_XICRA_PE(folder_path, reads_path, name, threads_given, debug_bool, database_folder, seqtk_bin, R2_reads):
    
    ## debugging messages
    if debug_bool:
        print ("\n********* XICRA PE analysis *********\n")
    
    ## send XICRA command
    ## create argparse with arguments provided to call XICRA prep
    output_folder_XICRA = os.path.join(folder_path, 'analysis')
    XICRA_options_prep = argparse.Namespace(input=reads_path, output_folder=output_folder_XICRA, 
                                            single_end=False, batch=False, in_sample=False, 
                                            ex_sample=False, detached=False, include_lane=False, 
                                            include_all=False, threads=threads_given, merge_Reads=False, 
                                            copy_reads=False, rename=False, help_format=False, 
                                            help_project=False, debug=debug_bool)
    prep.run_prep(XICRA_options_prep)
    
    ## create argparse with arguments provided to call XICRA prep
    XICRA_options_join = argparse.Namespace(input=output_folder_XICRA, noTrim=True,
                                            single_end=False, batch=False, in_sample=False, 
                                            ex_sample=False, detached=False, include_lane=False, 
                                            include_all=False, threads=threads_given, perc_diff= 8,  
                                            help_format=False,  help_project=False, help_join_reads=False, debug=debug_bool)
    join.run_join(XICRA_options_join)
    
    ## create argparse with arguments provided to call XICRA prep
    software_miRNA = ["sRNAbench", "optimir", "miraligner"]
    XICRA_options_miRNA = argparse.Namespace(input=output_folder_XICRA, 
                                            single_end=False, batch=False, in_sample=False, noTrim=False, 
                                            ex_sample=False, detached=False, include_lane=False, 
                                            include_all=False, threads=threads_given,
                                            soft_name=software_miRNA, species='hsa',
                                            database=database_folder, miRNA_gff=False, hairpinFasta=False, 
                                            matureFasta=False, miRBase_str=False, 
                                            help_format=False,  help_project=False, help_miRNA=False, debug=debug_bool)
    miRNA.run_miRNA(XICRA_options_miRNA)
    
    ## debugging messages
    if debug_bool:
        print ("\n********* XICRA R1 analysis *********\n")
    
    #########
    ## R1
    #########
    R1_in_file = os.path.join(folder_path, "R1.txt")
    with open(R1_in_file, 'w') as fh:
        fh.write("R1")
    fh.close()
    
    ## create argparse with arguments provided to call XICRA prep
    output_folder_XICRA_R1 = os.path.join(folder_path, 'analysis_R1')
    XICRA_options_miRNA_R1 = argparse.Namespace(input=reads_path, output_folder =output_folder_XICRA_R1,
                                            single_end=True, batch=False, in_sample=R1_in_file, 
                                            ex_sample=False, detached=True, include_lane=False, 
                                            include_all=False, threads=threads_given, noTrim=True,
                                            soft_name=software_miRNA, species='hsa',
                                            database=database_folder, miRNA_gff=False, hairpinFasta=False, 
                                            matureFasta=False, miRBase_str=False, 
                                            help_format=False,  help_project=False, help_miRNA=False, debug=debug_bool)
    miRNA.run_miRNA(XICRA_options_miRNA_R1)
    
    
    #########
    ## R2
    #########
    ## rev_comp reads
    ## debugging messages
    if debug_bool:
        print ("\n********* XICRA R2 analysis *********\n")
    
    R2_reads_revComp = R2_reads.split("R2.fq")[0] + "revComp.fq" 
    seqtk_cmd = "%s seq -r %s > %s" %(seqtk_bin, R2_reads, R2_reads_revComp)
    print ("+ Reverse complement reads")
    code = functions.system_call(seqtk_cmd)
    if not code:
        print ("** ERROR: Some error occurred when using seqtk for reverse complement")
        exit()
        
    ##
    R2_in_file = os.path.join(folder_path, "R2.txt")
    with open(R2_in_file, 'w') as fh2:
        fh2.write("revComp")
    fh2.close()
    
    output_folder_XICRA_R2 = os.path.join(folder_path, 'analysis_R2')
    XICRA_options_miRNA_R2 = argparse.Namespace(input=reads_path, output_folder =output_folder_XICRA_R2,
                                            single_end=True, batch=False, in_sample=R2_in_file, 
                                            ex_sample=False, detached=True, include_lane=False, 
                                            include_all=False, threads=threads_given, noTrim=True,
                                            soft_name=software_miRNA, species='hsa',
                                            database=database_folder, miRNA_gff=False, hairpinFasta=False, 
                                            matureFasta=False, miRBase_str=False, 
                                            help_format=False,  help_project=False, help_miRNA=False, debug=debug_bool)
    miRNA.run_miRNA(XICRA_options_miRNA_R2)
    return()

#########################################################
def process_fasta_length (fasta_file, folder, debug): 
    
    len_Dataframe = pd.DataFrame(columns=('len', 'ID' , 'seq'))
    
    ## read file    
    with open(fasta_file, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 2:
                record = reads2tabular.process_fasta(lines)
                # re-init
                lines = []
                len_Seq = len(record['sequence'])
                len_Dataframe.loc[ len(len_Dataframe) ] = (len_Seq, record['name'], record['sequence']) 

    ##
    grouped_df = len_Dataframe.groupby(['len'])
    len_dict = {}
    
    for len_int, cluster in grouped_df:
        ## write file    
        file_name = os.path.join(folder, 'seqs_len' + str(len_int) + '.fa')
        
        ## debugging messages
        if debug:
            print ("** Printing reads of length (%s) in file %s" %(len_int, file_name))
        
        with open(file_name, 'w') as outfh:
            for index, row in cluster.iterrows():
                outfh.write(row['ID'] + '\n' + row['seq'] + '\n')
        outfh.close()
        len_dict [len_int] = file_name
        
    return (len_dict)

#####################################################
parser = argparse.ArgumentParser(prog='simulation_sender.py', formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description='''

simulation_sender.py: Given a fasta file and simulation parameters it 
                      generates a NGS reads simulation. It also sends XICRA miRNA command

Version: 0.1
License: GPLv3

USAGE: python simulation_sender.py --folder out_name --fasta fasta_file --reads PE SE 
                   --seqSys HS10 GA2 --art_bin path/art_illumina -l 16 --fcov 10 20 30 -t 4
                   [--freqs table.freqs.csv] [--n_rows 100] [--replicates 100]
                   [--debug]
              

 NOTE: sequencing system ID names are:
------------------------------------------------------------------------------------------------------------------------------
 GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)
------------------------------------------------------------------------------------------------------------------------------

''', epilog="Original code: JFSanchezHerrero")

#####################################################
parser.add_argument('--fasta', action='store', help='Fasta file to simulate', required=True)
parser.add_argument('--folder', action='store', help='Folder to store results', required=True)
parser.add_argument('--reads', dest='type_reads', nargs='*', help="Type of reads to simulate", choices=['PE', 'SE'], required=True )
parser.add_argument('--seqSys', dest='seqSys_list', nargs='*',
                    choices=['GA1', 'GA2', 'HS10', 'HS20', 'HS25', 'HSXn', 'HSXt', 'HiSeqX', 'MinS', 'MSv1', 'MSv3', 'NS50'], 
                    help='The name of Illumina sequencing system of the built-in profile used for simulation', required=True)

parser.add_argument('--art_bin', action='store', help='ART NGS simulation binary file [art_illumina]', required=True)
parser.add_argument('--seqtk_bin', action='store', help='seqtk binary file', required=True)
parser.add_argument('--fcov', dest='fcov_list', nargs='*', help='Fold coverage')
parser.add_argument('-t', '--threads', type=int, default=2)

parser.add_argument('--freqs', action='store', help='Frequencies table to subset nrows')
parser.add_argument('-n', '--n_rows', type=int, default=10)
parser.add_argument('-r', '--replicates', type=int, default=10)
parser.add_argument('--database', action='store', help="XICRA miRNA database or new folder for downloading necessary files")


parser.add_argument('--send_XICRA', action='store_true', default=True, help='Send XICRA analysis for each replicate')
parser.add_argument('--debug', action='store_true', default=False, help='Developer messages')
args = parser.parse_args()
#####################################################

print ("\n\n#############################################################################")
print ("# Simulate NGS reads for sequence fasta file: " + args.fasta)
print ("# Type of reads to simulate:: " + str(args.type_reads))
print ("# Using built-in profile: " + str(args.seqSys_list))
print ("# Read length: Determined by isomiR length")
print ("# Fold coverage: " + str(args.fcov_list))
print ("# Output folder: " + os.path.abspath(args.folder))

if (args.freqs):
   print ("# Subset frequencies table provided: " + os.path.abspath(args.freqs))
   print ("# Number of rows to subset for replicate: " + str(args.n_rows))
   print ("# Number of replicates to do: " + str(args.replicates))

print ("#############################################################################\n")

## create main folder
args.folder = os.path.abspath(args.folder)
functions.create_folder(args.folder)

## main script path
dirnamePath=os.path.dirname(sys.argv[0])

## os.path.abspath
args.fasta = os.path.abspath(args.fasta)
args.art_bin = os.path.abspath(args.art_bin)
args.freqs = os.path.abspath(args.freqs)

if (args.freqs):
    ## create X subsets (replicates) using N (n_rows) miRNA each time
    ## simulate NGS reads for each replicate
    for i in range(args.replicates):
        str_rep = "rep_" + str(i + 1)
        print ("****************************************************************************")
        print ("+ Replicate: " + str_rep)
        replicate_path = functions.create_subfolder(str_rep, args.folder)
        
        filename_stamp = replicate_path + '/.success'
        if os.path.isfile(filename_stamp):
            stamp =    functions.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s]" %(stamp, str_rep), 'yellow'))
        else:
            ## subset freqs table
            mod_freq_script = os.path.join(dirnamePath, "mod_freq.py")
            mod_freqs_file = os.path.join(replicate_path, str_rep + ".freqs")
            mod_freq_python_cmd = "python %s --freq %s --out %s --random_rows %s" %(
                mod_freq_script, args.freqs, mod_freqs_file, args.n_rows)
            
            ## function system command
            print ("+ Create random subset of miRNA freqs")
            code= functions.system_call(mod_freq_python_cmd)
            if not (code):
                print ("** ERROR: Something happened and the script failed...")
                exit()
                
            ## subset isomiRs
            get_isomiRs_script = os.path.join(dirnamePath, "get_isomiRs.py")
            mod_freqs_file_isomiRs = mod_freqs_file + '.isomiRs'
            get_isomiRs_python_cmd = "python %s --freq %s --out %s --fasta %s" %(
                get_isomiRs_script, mod_freqs_file + '.csv', mod_freqs_file_isomiRs, args.fasta)
            
            print ("+ Select miRNA sequences")
            code2 = functions.system_call(get_isomiRs_python_cmd)
            if not (code2):
                print ("** ERROR: Something happened and the script failed...")
                exit()
            
            ## simulate
            print ("+ Simulate NGS reads from miRNA sequences")
            subset_fasta = mod_freqs_file_isomiRs + '.fasta'
            NGS_simulator(str_rep, replicate_path, args.seqSys_list, args.type_reads, 
                          args.fcov_list, subset_fasta, args.threads, args.debug, 
                          args.art_bin, args.seqtk_bin, args.database, args.send_XICRA)
            
            print ("\n\n")
else:
    NGS_simulator('sim', os.path.abspath(args.folder), args.seqSys_list, 
                  args.type_reads, args.fcov_list, args.fasta, args.threads, 
                  args.debug, args.art_bin, args.seqtk_bin, args.database, args.send_XICRA)
