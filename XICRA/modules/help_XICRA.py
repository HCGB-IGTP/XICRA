#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain      ##
##########################################################
"""
Help messages for different scripts, modules
"""

###############
def help_fastq_format():
    """
    Explanation of fastq format details.
    
    See additional information in section under user-guide: 
    :ref:`format fastq files<format-fastq-files>` 
    """

    functions.boxymcboxface("Name format for samples")

    print ("Format for fastq files can be:")
    print ("name.fastq.gz, name_1.fastq.gz, name_R2.fastq.gz, name_L001_R1.fastq.gz, etc.")
    print ("\nThere are many options and here we provide some guidelines of the format.")
    print ("\n")

    functions.print_sepLine("*",20,"red")
    print ("[1] Length limitation")
    functions.print_sepLine("*",20,"red")
    print ("There is a limitation for the sample ID ('name') of 10 characters.")
    print (colored("** BacterialTyper provides an option to rename samples if necessary: module prep option --rename **", 'yellow'))
    print ("\n")

    functions.print_sepLine("*",20,"red")
    print ("[2] Single end files")
    functions.print_sepLine("*",20,"red")
    print ("It is possible to provide NGS single-end files although some steps ")
    print ("of the process could not be accomplished using single-end files.")
    print (colored('** Use option --single-end in the different BacterialTyper modules **', 'yellow'))
    print ("name.fastq.gz")
    print ("name.fastq")
    print ("name.fq")
    print ("\n")

    functions.print_sepLine("*",20,"red")
    print ("[3] Paired-end files")
    functions.print_sepLine("*",20,"red")
    print ("Paired-end files are full supported. The format for these files are:")
    print ("Read1 => name_1.fastq.g or name_R1.fastq.gz")
    print ("Read2 => name_2.fastq.gz or name_R2.fastq.gz")
    print ("\n")

    functions.print_sepLine("*",55,"red")
    print ("[4] Lane information:")
    functions.print_sepLine("*",55,"red")
    print ("In some cases, files might contain lane information (*L00x* and/or *00x*).")
    print ("BacterialTyper supports these names as long as follow these examples:")
    print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
    print ("name_L00x_1.fastq.gz\tname_L00x_2.fastq.gz")
    print ("name_L00x_R1_00x.fastq.gz\tname_L00x_R2_00x.fastq.gz")
    print (colored("\n** If you need to merge fastq files from different lanes, use option --merge within module prep **", 'yellow'))
    print ("\n")

    functions.print_sepLine("*",15,"red")
    print ("[5] Extensions:")
    functions.print_sepLine("*",15,"red")
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
