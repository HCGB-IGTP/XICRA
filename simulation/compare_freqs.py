#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
Get frequence of reads for each type, variant, etc
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import pandas as pd
import argparse
import numpy as np
from termcolor import colored
from collections import defaultdict

## import my modules
from XICRA.scripts import functions
from XICRA.scripts import reads2tabular

#####################################
def observed_data_analysis(observed_freqs):
        
    print ("# Read observed raw data")
    observed_counts = functions.get_data(observed_freqs, ',', '')
    observed_counts[['miRNA','variant', 'UID']] = observed_counts.ID.str.split("&",expand=True,)
    
    ## debugging messages
    if args.debug:
        print ("\n***************** Debug **********************")
        print ("observed_counts")
        print (observed_counts)
    
    ## read sequence isomiRs identified
    print ("# Read sequence data observed")
    observed_seqs_file = observed_freqs.split('.csv')[0] + '_seq.csv'
    observed_seqs = functions.get_data(observed_seqs_file, ',', 'index_col=0')
    
    ## debugging messages
    if args.debug:
        print ("\n***************** Debug **********************")
        print ("observed_seqs")
        print (observed_seqs)
    
    return(observed_counts, observed_seqs)

###################
def count_miRNA_fastq (fastq_file):
    fastq_count = defaultdict(int)
    with open(fastq_file, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 4:
                record = reads2tabular.process_fastq(lines)
                #sys.stderr.write("Record: %s\n" % (str(record)))
                lines = []
                fastq_ID = record['name'].split('::')[0].replace('@', '')
                fastq_count[fastq_ID] += 1
                
    return fastq_count        

###################
def analysis_observed_expected(name, given_tag, counts_observed, count_R1_reads, expected_counts_given, seqs_observed, isomiRDict):

    results_df = pd.DataFrame(columns=("name", "miRNA", "variant", "sequence", "obs", "exp", "TP", "FP", "FN", "S", "P"))

    ## analysis
    total_count_observed = counts_observed[given_tag].sum()
    miRNA_grouped = counts_observed.groupby(['miRNA'])
    for miRNA_ID, cluster in miRNA_grouped:
        total_count_observed_cluster = cluster[given_tag].sum()    
        total_count_expected = count_R1_reads[miRNA_ID.lower()]
        if args.debug:
            print ("###################")
            print ("miRNA: " + miRNA_ID)
            print (cluster)
            print ("Total observed:" + str(total_count_observed_cluster))
            print ("Total expected:" + str(total_count_expected))
            print ("###################")
            
        for index, row in cluster.iterrows():
            miRNA_variant = row['variant'] 
            observed_count = row[given_tag]
            if (',' in str(miRNA_variant)):
                ## complex variant
                if args.debug:  
                    print ("+ ---------------------------------------------------- +")
                    print ("** Complex variant: " + miRNA_variant)
                    print ("False Positive") ## complex variant
                    print ("observed_count:" + str(observed_count))
                    
                    ## TODO: check if any sequence expected matches the UID observed
                    ## see as example: hsa-miR-3944-3p,TTCGGGCTGGCCTGCTGCTCCGGGG
                    expected_count = 0
                    (TP, FP, FN, S, P) = get_results(observed_count, expected_count)
                    ## Do not save this result

            else:
                if args.debug:
                    print ("\n+ ---------------------------------------------------- +")
                    print ("miRNA: " + miRNA_ID)
                    print ("Variant: " + miRNA_variant)
                
                ## skip
                if (miRNA_variant == 'notsure'):
                    continue
                
                if (':' in miRNA_variant): 
                    if ('add' in miRNA_variant): ## iso_add5p, iso_add3p
                        miRNA_variant = miRNA_variant.split(':')[0]
                    else:
                        miRNA_variant = miRNA_variant[:-1] ## iso_5p, iso_3p, deletion
                        
                if (variantDict[miRNA_variant]): ## exists in variantDict
                    
                    if (miRNA_ID.lower() not in expected_counts_given.index):
                        continue
                        ## Do not save this result
                        #(TP, FP, FN, S, P) = get_results(observed_count, expected_count)
                        
                    entry_thisVariant = expected_counts_given.loc[miRNA_ID.lower(), variantDict[miRNA_variant]]
                    
                    if type(entry_thisVariant) is float:
                        entry_thisVariant_count = 0
                        entry_thisVariant_miRNA = ""
                    else:
                        entry_thisVariant = entry_thisVariant.split('x')
                        entry_thisVariant_count = int(entry_thisVariant[-1])
                        entry_thisVariant_miRNA = entry_thisVariant[0]
                    
                    ## count
                    expected_count=int(entry_thisVariant_count/100*total_count_expected)
                    
                    ## same sequence??
                    if (entry_thisVariant_miRNA in isomiRDict.keys()):
                        expected_seq_isomiR = isomiRDict[entry_thisVariant_miRNA]['seq']
                    else:
                        expected_seq_isomiR = ""
                    ##
                    observed_seq_isomiR = seqs_observed.loc[row['UID'], 'Read']
    
                    ## debugging messages
                    if args.debug:  
                        print ("\nVariant conversion:")
                        print (variantDict[miRNA_variant])
                        print ("entry_thisVariant_freq_count: " + str(entry_thisVariant_count))
                        print ("entry_thisVariant_miRNA: " + entry_thisVariant_miRNA)
                        print ("\nCounts\nexpected_count: " + str(expected_count))
                        print ("observed_count:" + str(observed_count))
                        print ("\nSequences\nLicense plate: " + row['UID'])
                        print ("observed_seq_isomiR: " + observed_seq_isomiR)
                        print ("expected_seq_isomiR: " + expected_seq_isomiR)
                        
                    
                    if (expected_seq_isomiR == observed_seq_isomiR):
                        (TP, FP, FN, S, P) = get_results(observed_count, expected_count)
                    elif (expected_seq_isomiR == ""):
                        (TP, FP, FN, S, P) = get_results(observed_count, expected_count)
                    else:
                        if args.debug:  
                            print ("\n******** ATTENTION **********")
                            print ("Same variant different isomiR")
                            print ("********************************\n")
                            print ("+ ---------------------------------------------------- +\n")
                        continue
                
                
                    ### Save results            
                    if args.debug:  
                        print ("save results for this entry")
                    #columns=("name", "miRNA", 'variant', "sequence", "obs", "exp", "TP", "FP", "FN", "S", "P"))
                    results_df.loc[len(results_df)] = (name, miRNA_ID, miRNA_variant, 
                                                       observed_seq_isomiR, observed_count, expected_count,
                                                       TP, FP, FN, S, P)
                
                else:
                    ## debugging messages
                    if args.debug:  
                        print ("Unknown variant") ## complex variant
                        print ("observed_count:" + str(observed_count))
                    
                    ## Do not save this result
                    #(TP, FP, FN, S, P) = get_results(observed_count, expected_count)
            
            ### debugging messages
            if args.debug:  
                print ("+ ---------------------------------------------------- +\n")

            

    ##
    return(results_df)

def get_results(observed_count, expected_count):
    ##
    TP=0
    FP=0
    FN=0
    if (observed_count>expected_count):
        TP=expected_count
        FP=abs(observed_count-expected_count)
        FN=0
    else:
        TP=observed_count
        FP=0
        FN=abs(expected_count-observed_count)

    ## sensitivity & precision
    if (TP==0):
        S=0
        P=0
    else:
        S=TP/(TP+FN)
        P=TP/(TP+FP)
    
    if args.debug:
        print ("Statistics:")  
        print ("TP: " + str(TP))
        print ("FP: " + str(FP))
        print ("FN: " + str(FN))
        print ("S: " + str(S))
        print ("P: " + str(P))

    return(TP, FP, FN, S, P)

#####################################################
parser = argparse.ArgumentParser(prog='compare_freqs.py', formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description='''

compare_freqs.py: Given the original freqs table (expected) and the simulated freqs 
(observed), determines precision, recall etc.

Version: 0.1
License: GPLv3

''', epilog="Original code: JFSanchezHerrero")

#####################################################
parser.add_argument('--name', action='store', help='Name to output results', required=True)
parser.add_argument('--folder', action='store', help='Folder containing full analysis', required=True)
parser.add_argument('--tag', action='store', help='Name of the sample')
parser.add_argument('--replicates', action='store_true', help='Analysis contains multiple replicates')

subparser_observed_freqs_name = parser.add_argument_group("Observed frequencies")
subparser_observed_freqs = subparser_observed_freqs_name.add_mutually_exclusive_group(required= True)
subparser_observed_freqs.add_argument('--observed_freqs', action='store', 
                                      help='Observed frequencies table (XICRA results)')
subparser_observed_freqs.add_argument('--retrieve_all', action='store_true', 
                                      help='Use folder provided to retrieve all results available (sRNAbench, miraligner, optimiR & PE, R1 and R2)')

parser.add_argument('--debug', action='store_true', default=False, help='Developer messages')
args = parser.parse_args()
#####################################################

#===============================================================================
# ----------------------------------
# Simulated isomiR categories: isomiR-Benchmark ##
# ----------------------------------
# FA <- Five add
# FS <- Five del
# NT <- Non-template
# SR <- SNP Rest
# SS <- SNP Seed
# TA <- Three add
# TS <- Three del
# CN <- Canonical/mature
#
# ----------------------------------
# Variant types for miRTop: 
# ----------------------------------
# iso_5p/iso_3p:+/-N. 
#  (+) indicates the start is shifted to the right. 
#  (-) indicates the start is shifted to the left. 
# N the number of nucleotides of difference. 
# For instance, if the sequence starts 2 nts after the reference miRNA, the label will be: iso_5p:+2, but if it starts before, 
# the label will be iso_5p:-2.
# 
# iso_add3p:N. Number of non-template nucleotides added at 3p.
# iso_add5p:N. Number of non-template nucleotides added at 5p.
# 
# iso_snv_seed: when affected nucleotides are between [2-7].
# iso_snv_central_offset: when affected nucleotides is at position [8].
# iso_snv_central: when affected nucleotides are between [9-12].
# iso_snv_central_supp: when affected nucleotides are between [13-17].
# iso_snv: anything else.
#
# ----------------------------------
# Conversion variant types 
# ----------------------------------
#
#    isomiR-Benchmark              miRTop
# ----------------------------------------------
#    Five add (FA)                iso_5p:-1
#    Three add (TA)               iso_3p:+1
#    Five del (FS)                iso_5p:+1
#    Three del (TS)               iso_3p:-1
#    Non-template (NT)            iso_add5p / iso_add3p
#    SNP Rest (SR)                iso_snv_central & iso_snv_central_supp
#    SNP Seed (SS)                iso_snv_seed & iso_snv_central_offset
#
# ----------------------------------
# Descriptive statistics 
# ----------------------------------
#
# Sensitivity = TP / (TP+FN)
# Precision = TP / (TP+FP)
#
# TP (true positives) indicates the number of counts that are correctly identified.
# FN (false negatives) indicates the number of genuine counts that were not detected.
# FP (false positives) indicates the number of simulated counts incorrectly identified.
#===============================================================================

## Dictionary for variant conversion: isomiR-Benchmkar vs miRTop
variantDict = {
    'NA':'CN',
    
    ####
    'iso_5p:-': 'FA',
    'iso_3p:+': 'TA',
    ####
    'iso_5p:+': 'FS',
    'iso_3p:-': 'TS',
    ####
    'iso_add5p': 'NT',
    'iso_add3p': 'NT',
    ####
    'iso_snv_seed': 'SS',
    'iso_snv_central_offset': 'SS',
    'iso_snv': 'SS',
    ####
    'iso_snv_central': 'SR',
    'iso_snv_central_supp': 'SR'
    
}

## main folder provided
folder = os.path.abspath(args.folder)
folder_rep_list = ()
rep_ID=""
## check how many replicates
if args.replicates:
    ## multiple replicates subfolders containing results
    print ("+ Analysis for multiple replicates provided")
    args.retrieve_all = True
    folder_rep_list = [os.path.join(folder, o) for o in os.listdir(folder) 
                       if os.path.isdir(os.path.join(folder,o))]
    print(folder_rep_list)
    ####
    
else:
    ## analysis for just one replicate
    print ("+ Analysis for one folder provided")
    folder_rep_list = [folder]
    print(folder)
    
    if args.retrieve_all:
        ## all analysis will be retrieved
        rep_ID= args.tag
    else:
        ##
        if not args.tag:
            print (colored("** ERROR: No option --tag provided. **", 'red'))
            exit()
        ##
        if '_revComp' in args.tag:
            rep_ID = args.tag.split('_revComp')[0]
        elif '_R1' in args.tag:
            rep_ID = args.tag.split('_R1')[0]
        else:
            rep_ID = args.tag

####
results = pd.DataFrame()
for folder_rep in folder_rep_list:
    if args.replicates:
        rep_ID = os.path.basename(folder_rep)
    
    ##
    print ("\n+ Analysis for: " + rep_ID)
    
    ## original counts
    print ("# Read expected frequency table")
    isomiRs_freqs = os.path.join(folder_rep, rep_ID + '.freqs.isomiRs.csv')
    isomiRs_fasta = os.path.join(folder_rep, rep_ID + '.freqs.isomiRs.fasta')
    expected_counts = functions.get_data(isomiRs_freqs, ',', 'index_col=0')
    
    ## debugging messages
    if args.debug:
        print ("\n***************** Debug **********************")
        print ("expected_counts")
        print (expected_counts)
    
    isomiR_dict =  {}
    col_list = list(expected_counts) ## get columns
    for miRNA, row in expected_counts.iterrows():
        for col in col_list:
            if (type(row[col]) == float):
               continue 
            if (row[col] == 0):
                continue
            entry = row[col].split('x')
            isomiR_dict[entry[0]] = {'count':entry[1], 'seq': ""} 
    
    ## debugging messages
    if args.debug:
        print ("\n***************** Debug **********************")
        print ("isomiR_dict")
        print (isomiR_dict)
    
    print ("# Read sequence data expected")
    with open(isomiRs_fasta, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 2:
                record = reads2tabular.process_fasta(lines)
                #sys.stderr.write("Record: %s\n" % (str(record)))
                lines = []
                fasta_ID = record['name'].replace('>', '')
                fasta_ID = fasta_ID.split('_')[0]
                if (fasta_ID in isomiR_dict.keys()):
                    isomiR_dict[fasta_ID]['seq'] = record['sequence'] 
    
    ## debugging messages
    if args.debug:
        print ("\n***************** Debug **********************")
        print ("isomiR_dict")
        print (isomiR_dict)
    
    ## get total count of reads
    total_count_miRNA = defaultdict(int)
    reads_R1 = os.path.join(folder_rep, 'reads', rep_ID + '_R1.fq')
    reads_R2 = os.path.join(folder_rep, 'reads', rep_ID + '_R2.fq')
    reads_R1_count = count_miRNA_fastq(reads_R1)
    reads_R2_count = count_miRNA_fastq(reads_R2)
    
    ## todo: check if count matches both reads
    
    ## debugging messages
    if args.debug:
        print ("\n***************** Debug **********************")
        print ("reads_R1_count")
        print (reads_R1_count)
        print ("reads_R2_count")    
        print (reads_R2_count)    
    
    ## get observed data
    if (args.observed_freqs):
        print ("+ Single file provided with observed data. Retrieve expected vs. observed results.")
        observed_counts, observed_seqs = observed_data_analysis(args.observed_freqs)
        results = analysis_observed_expected(args.name, args.tag, observed_counts, 
                                             reads_R1_count, expected_counts, 
                                             observed_seqs, isomiR_dict)
        continue
        
    elif (args.retrieve_all):
        print ("+ Retrieve all available comparisons from folder provided. Retrieve expected vs. observed results for all at the same time.")
        
        observed_counts_dict = {}
        
        ## PE
        PE_analysis_folder=os.path.join(folder_rep, 'analysis', 'report', 'miRNA')
        if os.path.isdir(PE_analysis_folder):
            files_PE_analysis = functions.retrieve_matching_files(PE_analysis_folder, ".csv")
            files_PE_analysis = [s for s in files_PE_analysis if '_seq' not in s]
            files_PE_analysis = [s for s in files_PE_analysis if '_dup' not in s]
            for f in files_PE_analysis:
                soft=f.split('expression-')[1].split('.csv')[0]
                observed_counts_dict[soft] = {'PE': f}
                
            
        ## R1
        R1_analysis_folder=os.path.join(folder_rep, 'analysis_R1', 'report', 'miRNA')
        if os.path.isdir(R1_analysis_folder):
            files_R1_analysis = functions.retrieve_matching_files(R1_analysis_folder, ".csv")
            files_R1_analysis = [s for s in files_R1_analysis if '_seq' not in s]
            files_R1_analysis = [s for s in files_R1_analysis if '_dup' not in s]
            for f in files_R1_analysis:
                    soft=f.split('expression-')[1].split('.csv')[0]
                    observed_counts_dict[soft]['R1'] =  f
                
        
        ## R2
        R2_analysis_folder=os.path.join(folder_rep, 'analysis_R2', 'report', 'miRNA')
        if os.path.isdir(R2_analysis_folder):
            files_R2_analysis = functions.retrieve_matching_files(R2_analysis_folder, ".csv")
            files_R2_analysis = [s for s in files_R2_analysis if '_seq' not in s]
            files_R2_analysis = [s for s in files_R2_analysis if '_dup' not in s]
            for f in files_R2_analysis:
                    soft=f.split('expression-')[1].split('.csv')[0]
                    observed_counts_dict[soft]['R2'] =  f

        ## debugging messages
        if args.debug:
            print ("\n***************** Debug **********************")
            print ("observed_counts_dict")
            print (observed_counts_dict)
            print()
        
        for soft_name in observed_counts_dict.keys():
            for type_read in observed_counts_dict[soft_name]:
                if type_read == 'R1':
                    tag_given = rep_ID + '_R1'
                elif type_read == 'R2':
                    tag_given = rep_ID + '_revComp'
                else:
                    tag_given = rep_ID

                ## debugging messages
                if args.debug:
                    print ("\n***************** Debug **********************")
                    print ("tag_given: " + tag_given)                    
                    print ("rep_ID: " + rep_ID)                    
                    print ("soft_name: " + soft_name)                    
                    print ("type_read: " + type_read)                    
                    print ("File: "+ observed_counts_dict[soft_name][type_read])
            
                print ("\t+ Analysis for: "+  rep_ID + " :: " + soft_name + " :: " + type_read)                    
                observed_counts, observed_seqs = observed_data_analysis(observed_counts_dict[soft_name][type_read])
                results_tmp = analysis_observed_expected(rep_ID, tag_given, 
                                                         observed_counts, reads_R1_count, expected_counts, 
                                                         observed_seqs, isomiR_dict)
                ## add soft_name and type_read
                results_tmp['soft'] = soft_name
                results_tmp['type_read'] = type_read
                
                if results.empty:
                    results=results_tmp
                else:
                    ## concat for all results
                    results = pd.concat([results, results_tmp])
        
##
print ("\n\n + Save simulation results in file:")
name = args.name + "_XICRA.simulations.csv"
print (name)
results.to_csv(name )
exit()
            