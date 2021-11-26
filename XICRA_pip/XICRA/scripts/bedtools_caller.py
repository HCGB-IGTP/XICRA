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

from XICRA.config import set_config
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.time_functions as HCGB_time

#import pybedtools

#######################################################
def subtract_coordinates(file1, file2, path_given, name, options, debug):
    """
    Function to create a substraction call using bedtools
    """
        
    ## create a name
    bed_file = os.path.join(path_given, name + ".bed") 
    #bed_file_tmp = bed_file + '_tmp' 
    
    ## check if previously done
    filename_stamp = path_given + '/.' + name + '_subtract_success'
    if os.path.isfile(filename_stamp):
        if HCGB_files.is_non_zero_file(bed_file):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'subtract bed annotations'), 'yellow'))
            return (bed_file)
    
    ## subtract any overlap between feature taking into account strandness
    string_options = " -s -f 0.5 -A " 
    if options:
        string_options = string_options + options
    
    ## Create call for bedtools intersect
    bedtools_exe = set_config.get_exe("bedtools", debug)
    cmd_bedtools = "%s subtract -a %s -b %s %s > %s" %(bedtools_exe, file1, file2, string_options, bed_file) 
    
    bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)
    if not bed_code:
        print(colored("** ERROR: Something happen while calling bedtools subtract for job: " + name, "red"))
        exit()
        
    ## print time stamp
    HCGB_time.print_time_stamp(filename_stamp)

    return (bed_file)    

#######################################################
def intersect_coordinates(file1, file2, path_given, name, options, debug):
    """
    Function to create intersect call using bedtools
    """
        
    ## create a name
    bed_file = os.path.join(path_given, name + ".bed") 
    #bed_file_tmp = bed_file + '_tmp' 
    
    ## check if previously done
    filename_stamp = path_given + '/.' + name + '_intersect_success'
    if os.path.isfile(filename_stamp):
        if HCGB_files.is_non_zero_file(bed_file):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'intersect bed annotations'), 'yellow'))
            return (bed_file)
    
    ###
    string_options = " -wa -wb "
    if options:
        string_options = string_options + options
    
    ## Create call for bedtools intersect
    bedtools_exe = set_config.get_exe("bedtools", debug)
    cmd_bedtools = "%s intersect -a %s -b %s %s > %s" %(bedtools_exe, file1, file2, string_options, bed_file) 
    
    bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)

    if not bed_code:
        print(colored("** ERROR: Something happen while calling bedtools intersect for job: " + name, "red"))
        exit()
        
    ## print time stamp
    HCGB_time.print_time_stamp(filename_stamp)

    return (bed_file)    

#######################################################
def convert_bam2bed(sample, bam_file, path_given, pilfer=False, debug=False):
    """
    This functions calls bedtools to generate a conversion from BAM to BED format.
    
    It converts, sorts and collapse information retaining counts for each feature in bed format
    """

    bed_file = os.path.join(os.path.abspath(path_given), sample + ".bed") ## create a name
    bed_file_tmp = bed_file + '_tmp' 
    
    filename_stamp = path_given + '/.convert_bam2bed_success'
    if os.path.isfile(filename_stamp):
        if HCGB_files.is_non_zero_file(bed_file):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, sample, 'convert_bam2bed'), 'yellow'))
            return (bed_file)
    
    ## execute conversion and count reads mapping in exact coordinates
    ## bedtools bamtobed -i bam_file > bed_file 
    ## bedtools sort -chrThenSizeA -i bed_file > bed_sort.file
    ## bedtools groupby -i bed_sort.file -o count -g 1,2,3 -c 4 > counts.bed 
    
    if pilfer:
        ## Create call in two separate calls to reduce RAM requirement
        bedtools_exe = set_config.get_exe("bedtools", debug)
        cmd_bedtools = "%s bamtobed -i %s > %s" %(bedtools_exe, bam_file, bed_file)
        bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)
    
        if bed_code:
            ## print time stamp
            HCGB_time.print_time_stamp(filename_stamp)
        
            return (bed_file)
        else:
            print ("** ERROR: Some error occurred during conversion from BAM to BED... **")
            exit()
    else:
        ## Create call in two separate calls to reduce RAM requirement
        bedtools_exe = set_config.get_exe("bedtools", debug)
        cmd_bedtools = "%s bamtobed -i %s | %s groupby -o count -g 1,2,3 -c 4 > %s" %(bedtools_exe, bam_file, bedtools_exe, bed_file_tmp)
        bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)
    
        if bed_code:
            cmd_bedtools2 = "%s sort -chrThenSizeA -i %s | %s groupby -o count -g 1,2,3 -c 4 > %s" %(bedtools_exe, bed_file_tmp, 
                                                                                                                    bedtools_exe, bed_file)
            bed_code2 = HCGB_sys.system_call(cmd_bedtools2, False, True)
        
        ## -----------------------------------------------
        ## Pybedtools
        ## -----------------------------------------------
        ## It might be possible to use pybedtools. We need to load bam, bed or whatever file. it is not possible to use string to absolute path
        #bed_info = pybedtools.bedtool.BedTool.bam_to_bed(bam_file)
        
        ## We might try but I guess to many RAM would be required.
        
        if not bed_code or not bed_code2:
            print ("** ERROR: Some error occurred during conversion from BAM to BED... **")
            exit()
    
    ## print time stamp
    HCGB_time.print_time_stamp(filename_stamp)

    ## remove tmp files
    os.remove(bed_file_tmp)

    return (bed_file)

## -----------------------------------------------------------
## Note on options of bedtools intersect
## -----------------------------------------------------------
# # Tool:    bedtools intersect (aka intersectBed)
# # Version: v2.29.2
# # Summary: Report overlaps between two feature files.
# # 
# # Usage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>
# # 
# #     Note: -b may be followed with multiple databases and/or 
# #     wildcard (*) character(s). 
# # Options: 
# #     -wa    Write the original entry in A for each overlap.
# # 
# #     -wb    Write the original entry in B for each overlap.
# #         - Useful for knowing _what_ A overlaps. Restricted by -f and -r.
# # 
# #     -loj    Perform a "left outer join". That is, for each feature in A
# #         report each overlap with B.  If no overlaps are found, 
# #         report a NULL feature for B.
# # 
# #     -wo    Write the original A and B entries plus the number of base
# #         pairs of overlap between the two features.
# #         - Overlaps restricted by -f and -r.
# #           Only A features with overlap are reported.
# # 
# #     -wao    Write the original A and B entries plus the number of base
# #         pairs of overlap between the two features.
# #         - Overlapping features restricted by -f and -r.
# #           However, A features w/o overlap are also reported
# #           with a NULL B feature and overlap = 0.
# # 
# #     -u    Write the original A entry _once_ if _any_ overlaps found in B.
# #         - In other words, just report the fact >=1 hit was found.
# #         - Overlaps restricted by -f and -r.
# # 
# #     -c    For each entry in A, report the number of overlaps with B.
# #         - Reports 0 for A entries that have no overlap with B.
# #         - Overlaps restricted by -f, -F, -r, and -s.
# # 
# #     -C    For each entry in A, separately report the number of
# #         - overlaps with each B file on a distinct line.
# #         - Reports 0 for A entries that have no overlap with B.
# #         - Overlaps restricted by -f, -F, -r, and -s.
# # 
# #     -v    Only report those entries in A that have _no overlaps_ with B.
# #         - Similar to "grep -v" (an homage).
# # 
# #     -ubam    Write uncompressed BAM output. Default writes compressed BAM.
# # 
# #     -s    Require same strandedness.  That is, only report hits in B
# #         that overlap A on the _same_ strand.
# #         - By default, overlaps are reported without respect to strand.
# # 
# #     -S    Require different strandedness.  That is, only report hits in B
# #         that overlap A on the _opposite_ strand.
# #         - By default, overlaps are reported without respect to strand.
# # 
# #     -f    Minimum overlap required as a fraction of A.
# #         - Default is 1E-9 (i.e., 1bp).
# #         - FLOAT (e.g. 0.50)
# # 
# #     -F    Minimum overlap required as a fraction of B.
# #         - Default is 1E-9 (i.e., 1bp).
# #         - FLOAT (e.g. 0.50)
# # 
# #     -r    Require that the fraction overlap be reciprocal for A AND B.
# #         - In other words, if -f is 0.90 and -r is used, this requires
# #           that B overlap 90% of A and A _also_ overlaps 90% of B.
# # 
# #     -e    Require that the minimum fraction be satisfied for A OR B.
# #         - In other words, if -e is used with -f 0.90 and -F 0.10 this requires
# #           that either 90% of A is covered OR 10% of  B is covered.
# #           Without -e, both fractions would have to be satisfied.
# # 
# #     -split    Treat "split" BAM or BED12 entries as distinct BED intervals.
# # 
# #     -g    Provide a genome file to enforce consistent chromosome sort order
# #         across input files. Only applies when used with -sorted option.
# # 
# #     -nonamecheck    For sorted data, don't throw an error if the file has different naming conventions
# #             for the same chromosome. ex. "chr1" vs "chr01".
# # 
# #     -sorted    Use the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input.
# # 
# #     -names    When using multiple databases, provide an alias for each that
# #         will appear instead of a fileId when also printing the DB record.
# # 
# #     -filenames    When using multiple databases, show each complete filename
# #             instead of a fileId when also printing the DB record.
# # 
# #     -sortout    When using multiple databases, sort the output DB hits
# #             for each record.
# # 
# #     -bed    If using BAM input, write output as BED.
# # 
# #     -header    Print the header from the A file prior to results.
# # 
# #     -nobuf    Disable buffered output. Using this option will cause each line
# #         of output to be printed as it is generated, rather than saved
# #         in a buffer. This will make printing large output files 
# #         noticeably slower, but can be useful in conjunction with
# #         other software tools and scripts that need to process one
# #         line of bedtools output at a time.
# # 
# #     -iobuf    Specify amount of memory to use for input buffer.
# #         Takes an integer argument. Optional suffixes K/M/G supported.
# #         Note: currently has no effect with compressed files.
# # 
# # Notes: 
# #     (1) When a BAM file is used for the A file, the alignment is retained if overlaps exist,
# #     and excluded if an overlap cannot be found.  If multiple overlaps exist, they are not
# #     reported, as we are only testing for one or more overlaps.


