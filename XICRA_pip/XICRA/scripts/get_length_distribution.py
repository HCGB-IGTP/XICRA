#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess
import argparse
from termcolor import colored
import pandas as pd

from XICRA.config import set_config
from XICRA.scripts import bedtools_caller

import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.time_functions as HCGB_time
import HCGB.format_conversion

#######################################################
def parse_bed_length(bed_file, sample, debug):

    #read a file
    fileReader = open(bed_file)
    
    header_list = [ "transcript", "biotype", "gene", "length", "chr", "count", "sample", "strand"]    
    print("\t".join(header_list))

    ## skip comments at the beginning of files
    while True:
        line = fileReader.readline()
        line = line.rstrip()
        field=line.strip().split('\t')

        ## example
        # 1    1218699    1218706    1    1    1217688    1218768    ENST00000494748    100    -    1217688    1218768    255,0,0    3    0,311,1837    2583,176,757    ENSG00000078808    retained_intron

        if line:
            
            ## ATTENTION: It takes many time and RAM to store results.
            ## Dump into file directly
            ## save data
            #data_parse.loc[field[7], 'sample'] = sample
            #data_parse.loc[field[7],'biotype'] = field[17]
            #data_parse.loc[field[7],'gene'] = field[16]
            #data_parse.loc[field[7],'length'] = abs(int(field[2])-int(field[1]))
            #data_parse.loc[field[7],'chr'] = field[0]
            #data_parse.loc[field[7],'count'] = field[3]
            #data_parse.loc[field[7],'strand'] = field[9]
            
            length_match = abs(int(field[2])-int(field[1]))
            
            ##          transcript  biotype    gene        length                           chr       count     sample  strand    
            list2print = [ field[7], field[17], field[16], str(length_match), str(field[0]), str(field[3]), sample, field[9] ]
            
            print("\t".join(list2print))
        else:
            fileReader.close()
            break

#######################################################
def get_length_dist(bam_file, GTF_info, name, chr, path_given, folder_GTF, debug=False):
    """
    Gets read length distribution from mapping file
    
    It initially converts BAM -> BED format, then checks if annotation file provided is splitted 
    and converted to BED format.
    
    It then intersects annotations provided and mapping reads and generates a table with statistics.

    :param bam_file: Mapping file in BAM format
    :param GTF_info: GTF annotation file
    :param name: Name of the sample
    :param chr: Option that states if Chr splitted is already done and GTF_info is a subset file
    :param path_given: Path to store results
    :param folder_GTF: Path to store splitted GTF files
    :param debug: Wether to show debug messages or not.
    
    :type bam_file: string
    :type GTF_info: string
    :type name: string
    :type chr: boolean
    :type path_given: string
    :type folder_GTF: string
    :type debug: boolean
    """
    
    ## debugging messages
    if debug:
        HCGB_aes.debug_message("******************")
        HCGB_aes.debug_message("get_length_dist function call: ", color="yellow")
        print()
        HCGB_aes.debug_message("bam_file: " + bam_file, color="yellow")
        HCGB_aes.debug_message("GTF_info: " + GTF_info, color="yellow")
        HCGB_aes.debug_message("name: " + name, color="yellow")
        HCGB_aes.debug_message("path_given: " + path_given, color="yellow")
        HCGB_aes.debug_message("folder_GTF: " + folder_GTF, color="yellow")
        HCGB_aes.debug_message("chr: " + str(chr), color="yellow")
        
    
    print("+ Converting BAM file to BED...")    
    
    ## converts BAM file in BED format
    mapping_bed_file = bedtools_caller.convert_bam2bed(name, bam_file, path_given, debug=debug)
    
    ## convert annotation in GTF to BED
    print("+ Converting GTF annotation file into BED...")
    if chr: 
        ## GTF_info is already splitted:
        (bed_files, gtf_files) = convert_GTF2bed(GTF_info, folder_GTF, num_files=0, debug=debug)
    else:    
        (bed_files, gtf_files) = convert_GTF2bed(GTF_info, folder_GTF, debug=debug)

    ## debugging messages
    if debug:
        HCGB_aes.debug_message("bed_files", color="yellow")
        print(bed_files)        
        print()
        HCGB_aes.debug_message("gtf_files", color="yellow")
        print(gtf_files)
        
    ## Loop for each reference sequence BED file
    ## and intersect with mapping BAM->BED file
    
    print("+ Intersecting GTF annotation file and mapping BED file...")

    original_stdout = sys.stdout # Save a reference to the original standard output

    for RefSeq_bed_file in bed_files.values():
        ## intersect
        ## call bedtools_caller.intersect_coordinates()
        bed_file = bedtools_caller.intersect_coordinates(mapping_bed_file, RefSeq_bed_file, path_given, name, "-s", debug)
        out_file = os.path.join(path_given, name + ".length.results.txt") 
    
        ## print results to file
        with open(out_file, 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            
            ## print results directly to open file
            parse_bed_length(bed_file, name, debug)
            f.close()
    
    ## Reset print standard out
    sys.stdout = original_stdout
    
    ## Maybe merge all data from different chromosomes
    
    ## call R script in XICRA.stats to plot results
    
    exit()
    
#######################################################
def convert_GTF2bed(GTF_file, path_given, cpus2use=2, num_files=1, debug=False):
    """
    This functions calls gtf2bed from HCGB to generate a conversion from GTF to BED format.
    
    It allows to split big GTF file provided into multiple files and convert each subfile into BED format.
    
    """
    
    HCGB_files.create_folder(path_given)
    
    split_GTF_files = {}
    if num_files==0:
        ## GTF file is already splitted
        split_GTF_files['GTF_file'] = os.path.abspath(GTF_file)
    else:
        
        print("+ Split GTF in multiple subsets to speed up process")
                
        ## split GTF file
        split_GTF_files = HCGB.format_conversion.file_splitter.split_file_call(GTF_file, num_files=1, 
                                                 name=False, chr_option=True, 
                                                 in_format="GTF",
                                                 path_given=path_given, debug=debug) 
    # create bed
    print("+ Converting GTF files in BED format:")
    bed_files = {}
    for each_file in split_GTF_files.values():
        ## create name
        each_file_name = HCGB_files.get_file_name(each_file)
        bed_file = os.path.join(path_given, each_file_name + ".bed") ## create a name
        
        print("\t- Converting GTF file: " + each_file)
        print("\t- To BED file: " + each_file)
        print()
        ## convert GTF
        HCGB.format_conversion.gtf2bed.parse_GTF_call(each_file, bed_file, debug)
        
        ## save
        bed_files[each_file] = bed_file
    
    return (bed_files, split_GTF_files)

#######################################################
def main():
    ## this code runs when call as a single script
    parser=argparse.ArgumentParser(description='''Get length distribution for biotypes mapped in BAM file
    Note:
    It takes into account no genes or transcript are broken
    It is also possible to split according to chromosome (one gtf/chromosome)
    
    ''');
    
    parser.add_argument('--input', '-i', help='Input BAM file', required=True);

    parser.add_argument('--name', '-n',
                        help='Name of the sample. Default: use filename provided.', default="");
    
    parser.add_argument('--path', '-p',
                        help='Path to save results generated. Default: use path from filename provided.', default="");

    parser.add_argument('--path_GTF',
                        help='Path to store results generated from split GTF file. Default: use path provided previously.', default="");

    parser.add_argument('--annot', '-a',
                        help='Absolute path for GTF file.', default="");

    parser.add_argument('--annot_splitted', '-as', action="store_true",
                        help='--annot is a dir containing splitted GTF files for each chromosome.');
    
    args=parser.parse_args();
    
    ## lets split the big file provided
    print()
    print("** get_length_distribution script call **")
    print()
    files_generated = get_length_dist(os.path.abspath(args.input), os.path.abspath(args.annot), name=args.name, 
              chr=args.annot_splitted, path_given=os.path.abspath(args.path), folder_GTF=os.path.abspath(args.path_GTF), 
              debug=False)
    return ()

######
if __name__== "__main__":
    main()
    