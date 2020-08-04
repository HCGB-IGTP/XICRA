#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain      ##
##########################################################
"""
Help messages for different scripts, modules
"""
from termcolor import colored

from HCGB import functions


###############
def help_fastq_format():
    """
    Explanation of fastq format details.
    """

    functions.aesthetics_functions.boxymcboxface("Name format for samples")

    print ("Format for fastq files can be:")
    print ("name.fastq.gz, name_1.fastq.gz, name_R2.fastq.gz, name_L001_R1.fastq.gz, name_L001_R1_001.fastq.gz etc.")
    print ("\nThere are many options and here we provide some guidelines on the name format.")
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("[1] Length limitation")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("There is a limitation for the sample ID ('name') of 25 characters.")
    print (colored("** XICRA provides an option to rename samples if necessary: module prep option --rename **", 'yellow'))
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("[2] Single end files")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print (colored('** Use option --single-end in the different XICRA modules **', 'yellow'))
    print ("name.fastq.gz")
    print ("name.fastq")
    print ("name.fq")
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("[3] Paired-end files")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Paired-end files are full supported. The format for these files are:")
    print ("Read1 => name_1.fastq.g or name_R1.fastq.gz")
    print ("Read2 => name_2.fastq.gz or name_R2.fastq.gz")
    print (colored('** See additional details for Lane information **', 'yellow'))
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("[4] Lane information:")
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("In some cases, files might contain lane information (*L00x* and/or *00x*).")
    print ("XICRA supports these names as long as follow these examples:")
    print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
    print ("name_L00x_1.fastq.gz\tname_L00x_2.fastq.gz")
    print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
    print ("name_L00x_R1_00x.fastq.gz\tname_L00x_R2_00x.fastq.gz")
    print ("\n")
    print ("Sometimes it might be appropriate to include lane tags (*L00X*) within the name.")
    print (colored("** Use option --include-lane within each module", 'yellow'))
    
    print (colored("\n** If you need to merge fastq files from different lanes, use option within module prep **", 'yellow'))
    print("As an example:")
    print (colored("\n** Option --merge within module prep **", 'yellow'))
    print ("sample1_L001_R1.fastq.gz\tsample1_L001_R2.fastq.gz")
    print ("sample1_L002_R1.fastq.gz\tsample1_L002_R2.fastq.gz")
    print ("sample1_L003_R1.fastq.gz\tsample1_L003_R2.fastq.gz")
    print ("sample1_L004_R1.fastq.gz\tsample1_L004_R2.fastq.gz")
    print ("Result:")
    print ("--------------------------------------------------")
    print ("sample1_R1.fastq.gz\tsample1_R2.fastq.gz")
    print ("\n")
    print (colored("\n** Option --merge-by-lane within module prep **", 'yellow'))
    print ("sample1_L001_R1_001.fastq.gz\tsample1_L001_R2_001.fastq.gz")
    print ("sample1_L001_R1_002.fastq.gz\tsample1_L001_R2_002.fastq.gz")
    print ("sample1_L002_R1_001.fastq.gz\tsample1_L002_R2_001.fastq.gz")
    print ("sample1_L002_R1_002.fastq.gz\tsample1_L002_R2_002.fastq.gz")
    print ("--------------------------------------------------")
    print ("Result:")
    print ("sample1_L001_R1.fastq.gz\tsample1_L001_R2.fastq.gz")
    print ("sample1_L002_R1.fastq.gz\tsample1_L002_R2.fastq.gz")
    print (colored("** Remember to use option --include_lane within each module", 'yellow'))
    print ("\n")
    
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("[5] Include all information:")
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("In some cases, files might contain other information and it is necessay to " +
           "include it all as a tag nane. See as an example:")
    print ("sample1_L001_XYZ_R1_001.fastq.gz\tsample1_L001_XYZ_R2_001.fastq.gz")
    print (colored("** Remember to use option --include_all within each module", 'yellow'))
    print (colored("** It might be appropiate to change samples names using --rename option under prep module", 'yellow'))
    
    print ("\n")
    functions.aesthetics_functions.print_sepLine("*",15,"red")
    print ("[6] Extensions:")
    functions.aesthetics_functions.print_sepLine("*",15,"red")
    print ("name_L00x_R2.fastq\tname_L00x_R2.fq\nname_L00x_R2.fastq.gz\tname_L00x_R2.fq.gz")
    print ("\n")


def project_help():
    return ()

def multiqc_help():
    return ()

def print_help_adapters():
    return()

def help_join_reads():
    return ()

def help_miRNA():
    return ()
