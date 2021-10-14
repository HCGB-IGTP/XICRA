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
    print ("name.fastq.gz")
    print ("name_1.fastq.gz, adding '1' or '2' to specify the read")
    print ("name_R2.fastq.gz, adding 'R1' or 'R2' to specify the read")
    print ("name_L001_R1.fastq.gz, adding the lane information as L00x after the name")
    print ("name_L001_R1_001.fastq.gz, adding 00X at the end. This naming is useful when the fastq")
    print ("files of the sample sample had been cut in different files).")
    print ("name_L001_XYZ_R1_001.fastq.gz, there can be extra info for each file.")
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
    print ("[4] Sample identification:")
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    
    print ("XICRA will store the names of all the input files. After that, it will identify the samples.")
    print ("It can be the case that more than one file belong to the same sample. In order to pass this information")
    print ("to XICRA a combination of the following parameters may be needed:")
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("[4.1] Lane information:")
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("In some cases, files might contain lane information (*L00x* and/or *00x*).")
    print ("XICRA supports these names as long as follow these examples:")
    print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
    print ("name_L00x_1.fastq.gz\tname_L00x_2.fastq.gz")
    #print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
    print ("name_L00x_R1_00x.fastq.gz\tname_L00x_R2_00x.fastq.gz")
    print ("\n")
    
    print ("- If you want to include lane tags (*L00X*) into each  each sample name (differentiate samples considering the lane):")
    print (colored("** Use option --include-lane within each module and the lane tag will also be used to identify samples", 'yellow'))
    
    print ("\n- However, if you want to consider as a single sample the different lanes, you need to merge")
    print ("the fastq files from the different lanes, use option --merge_Reads within module prep.")
    print("As an example:")
    print (colored("** Options --merge_Reads within module prep **", 'yellow'))
    print ("sample1_L001_R1.fastq.gz\tsample1_L001_R2.fastq.gz")
    print ("sample1_L002_R1.fastq.gz\tsample1_L002_R2.fastq.gz")
    print ("sample1_L003_R1.fastq.gz\tsample1_L003_R2.fastq.gz")
    print ("sample1_L004_R1.fastq.gz\tsample1_L004_R2.fastq.gz")
    print ("Result:")
    print ("--------------------------------------------------")
    print ("sample1_R1.fastq.gz\tsample1_R2.fastq.gz")
    print ("\n")
    
    print ("\n- If you need to merge fastq files of the same lane that differ in the last group of numbers")
    print ("use option --mergeReads together with --include-lane within module prep.")
    print (colored("\n** Option --include_lane --merge-by-lane within module prep **", 'yellow'))
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
    
    ### if you want to merge lane and extension --mergeReads
    print ("- If you need to merge fastq files with different lanes and final extension ")
    print ("(_001, _002, ...), use only option --merge_Reads within module prep.")
    print("As an example:")
    print (colored("\n** Options --merge_Reads within module prep **", 'yellow'))
    print ("sample1_L001_R1_001.fastq.gz\tsample1_L001_R2_001.fastq.gz")
    print ("sample1_L001_R1_002.fastq.gz\tsample1_L001_R2_002.fastq.gz")
    print ("sample1_L002_R1_001.fastq.gz\tsample1_L002_R2_001.fastq.gz")
    print ("sample1_L002_R1_002.fastq.gz\tsample1_L002_R2_002.fastq.gz")
    print ("--------------------------------------------------")
    print ("Result:")
    print ("sample1_R1.fastq.gz\tsample1_R2.fastq.gz")
    print ("\n")    
    
    
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("[4.2] Include all information:")
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("In some cases, files might contain other extra information and it is necessary to ")
    print ("include it all as a tag name, in that case use --include-all. In the following example")
    print ("XYZ is the extra information and it is also used to identify each sample:")
    print ("sample1_L001_XYZ_R1_001.fastq.gz\tsample1_L001_XYZ_R2_001.fastq.gz")
    print (colored("** Remember to use option --include_all within each module", 'yellow'))
    
    print (colored("** It might be appropriate to change samples names using --rename option under prep module", 'yellow'))
    
    print ("\n- If you need to merge fastq files that only differ in the last group of numbers ")
    print ("(_001, _002, ...), use option --merge_Reads within module prep together with --include-all.")
    print("As an example:")
    print (colored("\n** Options --include_all --merge_Reads within module prep **", 'yellow'))
    print ("sample1_L001_XYZ_R1_001.fastq.gz\tsample1_L001_XYZ_R2_001.fastq.gz")
    print ("sample1_L001_XYZ_R1_002.fastq.gz\tsample1_L001_XYZ_R2_002.fastq.gz")
    print ("sample1_L002_XYZ_R1_001.fastq.gz\tsample1_L002_XYZ_R2_001.fastq.gz")
    print ("sample1_L002_XYZ_R1_002.fastq.gz\tsample1_L002_XYZ_R2_002.fastq.gz")
    print ("--------------------------------------------------")
    print ("Result:")
    print ("sample1_L001_XYZ_R1.fastq.gz\tsample1_L001_XYZ_R2.fastq.gz")
    print ("sample1_L002_XYZ_R1.fastq.gz\tsample1_L002_XYZ_R2.fastq.gz")
    print (colored("** Remember to use option --include_all within each module", 'yellow'))
    print ("\n")
    
    print ("\n")
    functions.aesthetics_functions.print_sepLine("*",15,"red")
    print ("[4.3] Extensions:")
    functions.aesthetics_functions.print_sepLine("*",15,"red")
    print ("name_L00x_R2.fastq\tname_L00x_R2.fq\nname_L00x_R2.fastq.gz\tname_L00x_R2.fq.gz")
    print ("\n")

###############
def project_help():
    return ()

###############
def multiqc_help():
    return ()

###############
def print_help_adapters():
    return()

###############
def help_join_reads():
    return ()

###############
def help_miRNA():
    return ()

###############
def help_tRNA():
    return ()

