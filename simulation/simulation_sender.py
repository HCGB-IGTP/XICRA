#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
from builtins import len
from numpy.f2py.auxfuncs import debugcapi
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

## import my modules
from XICRA.scripts import functions
from XICRA.scripts import reads2tabular

## paired - end y luego hacer R1 y R2
## HS25 HiSeq 2500
## cobertura: 10

#########################################################
def NGS_simulator(name, abs_folder, seqSys_list, type_reads, fcov_list, fasta, threads, debug, art_illumina_bin):

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
                reads_path = abs_folder
            else:
                reads_path = functions.create_subfolder(reads, profile_path)

            for fcov in fcov_list:
                print ("   + Coverage x" + fcov)

                ## check how many profiles to use
                if len(type_reads) == 1:
                    coverage_path = abs_folder
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
                    
                    ## add random seed? --rndSeed
                    
                    ## outfile
                    outfile_name = name + '_P-' + profile + '_T-' + reads + '_x-' + fcov + '_L-' + str(int_len) + '_R'
                    outfile_path = os.path.join(tmp_fastq, outfile_name)
                    
                    ## command
                    art_illumina_cmd = "%s -na -p -m 50 -s 5 -ss %s -i %s -l %s -f %s -o %s" %(art_illumina_bin, profile, fasta_file_len, str_len, fcov, outfile_path)
                    
                    code = functions.system_call(art_illumina_cmd)
                    if not code:
                        print ("** ERROR: Some error happened during ART simulation")
                        exit()
                        
                ## send XICRA command

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

USAGE: python simulation_sender.py --foolder out_name --fasta fasta_file --reads PE SE 
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
parser.add_argument('--fcov', dest='fcov_list', nargs='*', help='Fold coverage')
parser.add_argument('-t', '--threads', type=int, default=2)

parser.add_argument('--freqs', action='store', help='Frequencies table to subset nrows')
parser.add_argument('-n', '--n_rows', type=int, default=10)
parser.add_argument('-r', '--replicates', type=int, default=10)

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
        
        ## subset freqs table
        mod_freq_script = os.path.join(dirnamePath, "mod_freq.py")
        mod_freqs_file = os.path.join(replicate_path, str_rep + ".freqs")
        mod_freq_python_cmd = "python %s --freq %s --out %s --random_rows %s" %(mod_freq_script, args.freqs, mod_freqs_file, args.n_rows)
        ## function system command
        print ("+ Create random subset of miRNA freqs")
        code= functions.system_call(mod_freq_python_cmd)
        if not (code):
            print ("** ERROR: Somethin happened and the script failed...")
            exit()
            
        ## subset isomiRs
        get_isomiRs_script = os.path.join(dirnamePath, "get_isomiRs.py")
        mod_freqs_file_isomiRs = mod_freqs_file + '.isomiRs'
        get_isomiRs_python_cmd = "python %s --freq %s --out %s --fasta %s" %(get_isomiRs_script, mod_freqs_file + '.csv', mod_freqs_file_isomiRs, args.fasta)
        print ("+ Select miRNA sequences")
        code2 = functions.system_call(get_isomiRs_python_cmd)
        if not (code2):
            print ("** ERROR: Somethin happened and the script failed...")
            exit()
        
        ## simulate
        print ("+ Simulate NGS reads from miRNA sequences")
        subset_fasta = mod_freqs_file_isomiRs + '.fasta'
        NGS_simulator(str_rep, replicate_path, args.seqSys_list, args.type_reads, args.fcov_list, subset_fasta, args.threads, args.debug, args.art_bin)
        print ("\n\n")
else:
    NGS_simulator('sim', os.path.abspath(args.folder), args.seqSys_list, args.type_reads, args.fcov_list, args.fasta, args.threads, args.debug, args.art_bin)
