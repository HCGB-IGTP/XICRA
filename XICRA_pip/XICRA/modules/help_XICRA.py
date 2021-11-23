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


    print ("Format for fastq files can be:\n")
    print ("- name.fastq.gz")
    print ("- name_1.fastq.gz, with '1' or '2' to specify the read")
    print ("- name_R2.fastq.gz, 'R1' or 'R2' to specify the read")
    print ("- name_L001_R1.fastq.gz, with the lane information as 'L00X' or '00X' after the name")
    print ("- name_L001_R1_001.fastq.gz, with '00X' at the end of the file name. This naming is used") 
    print ("- when the fastq of a sample had been cut in different files.")
    print ("- name_L001_XYZ_R1_001.fastq.gz, there can be extra info for each file, XYZ.")
    print ("\nThe input file names should be structured considering the following aspects:")
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Length limitation")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("There is a limitation for the sample ID ('name') of 25 characters.")
    print (colored("** XICRA provides an option to rename samples if necessary: module prep option --rename **", 'yellow'))
    print ("\n")
    
    functions.aesthetics_functions.print_sepLine("*",15,"red")
    print ("Extensions:")
    functions.aesthetics_functions.print_sepLine("*",15,"red")
    print("The suported extensions are:\n")
    print ("- name_L00x_R2.fastq\tname_L00x_R2.fq\n- name_L00x_R2.fastq.gz\tname_L00x_R2.fq.gz")
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Single-end files")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print("It is possible to provide NGS single-end files although some steps of the process could not be accomplished")
    print("using single-end files.\n")    
    print ("- name.fastq.gz")
    print ("- name.fastq")
    print ("- name.fq")
    print (colored('** Use option --single-end in the different XICRA modules. **', 'yellow'))
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Paired-end files")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Paired-end files are full supported. The format for these files are:\n")
    print ("- name_1.fastq.gz, name_2.fastq.gz")
    print ("- name_R1.fastq.gz, name_R2.fastq.gz")
    print (colored('** No parameter is needed in to specify this kind of files. **', 'yellow'))
    print ("\n")    
    
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Lane information")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print("Files might contain lane information (L00x and/or 00x). XICRA")
    print("supports these names as long as follow these examples:")
    print("- name_L00x_R1.fastq.gz, name_L00x_R2.fastq.gz")
    print("- name_L00x_1.fastq.gz, name_L00x_2.fastq.gz")
    print ("\n")
    
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Name extensions")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print("It can also be the case that the reads of a sample are divided in different files.")
    print("In those cases, the files should contain a name final extension: ")
    print("- name1_L001_R1_001.fastq.gz, name1_L001_R2_001.fastq.gz")
    print("- name1_L001_R1_002.fastq.gz, name1_L001_R2_002.fastq.gz")
    print("- name1_L002_R1_001.fastq.gz, name1_L002_R2_001.fastq.gz")
    print("- name1_L002_R1_002.fastq.gz, name1_L002_R2_002.fastq.gz")
    print ("\n")
    
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print ("Extra information")
    functions.aesthetics_functions.print_sepLine("*",20,"red")
    print("In some cases, files might contain other extra information. In the following example,")
    print("XYZ is the extra information:")
    print("- name1_L001_XYZ_R1_001.fastq.gz, name1_L001_XYZ_R2_001.fastq.gz")
    print("- name1_L001_XYZ_R1_002.fastq.gz, name1_L001_XYZ_R2_002.fastq.gz")
    print ("\n")
    
    functions.aesthetics_functions.boxymcboxface("Sample identification")
        
    print ("XICRA will store the names of all the input files. After that, it will identify the samples.")
    print ("It can be the case that more than one file belong to the same sample. In order to pass this information")
    print ("to XICRA, a combination of the following parameters may be needed depending on the characteristics of the")
    print ("input file names:")
    print ("\n")

    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("Option --include_lane:")
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("If you want to include lane tags (L00X, 00X) into each  each sample name (differentiate samples considering the lane):")
    print (colored('** Use option --include_lane within each module and the lane tag will also be used to identify samples. **\n', 'yellow'))
    print("However, if you want to consider as a single sample the different lanes, you need to merge the")
    print("corresponding fastq files:")
    print (colored('** Use option --merge_Reads within module prep. **\n', 'yellow'))
    
    print("As an example, considering the input files:")

    print("- name1_L001_R1.fastq.gz, name1_L001_R2.fastq.gz")
    print("- name1_L002_R1.fastq.gz, name1_L002_R2.fastq.gz")
    print("- name1_L003_R1.fastq.gz, name1_L003_R2.fastq.gz")
    print("- name1_L004_R1.fastq.gz, name1_L004_R2.fastq.gz\n")
    
    print("\t 1.By adding the option --include_lane in all modules, XICRA will identify four samples:")    
    print("\t - Sample 1: name1_L001_R1, name1_L001_R2")
    print("\t - Sample 2: name1_L002_R1, name1_L002_R2")
    print("\t - Sample 3: name1_L003_R1, name1_L003_R2")
    print("\t - Sample 4: name1_L004_R1, name1_L004_R2")
    print (colored('\t ** Remember to use option --include_lane within each module. **\n', 'yellow'))
    
    print("\t 2. By adding the options --include_lane --merge_Reads within module prep, XICRA will only")
    print("\t identify one sample, merging all the corresponding files:")
    print("\t - Sample 1: sample1_R1, sample1_R2\n")
    
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print ("Option --include_all:")
    functions.aesthetics_functions.print_sepLine("*",55,"red")
    print("In some cases, files might contain other extra information and it is necessary to use all")  
    print("file name to identify samples:")
    print (colored('** If that is the case use --include_all in al modules. **\n', 'yellow'))
    print("If you want to merge fastq files that only differ in the final extension (_001, _002, ...):")
    print (colored('** Use options --merge_Reads --include_all within module prep and only --include_all in the rest of the modules. **\n', 'yellow'))
    
    print("As an example, considering the input files:")
    print("- name1_L001_XYZ_R1_001.fastq.gz, name1_L001_XYZ_R2_001.fastq.gz")
    print("- name1_L001_XYZ_R1_002.fastq.gz, name1_L001_XYZ_R2_002.fastq.gz")
    print("- name1_L002_XYZ_R1_001.fastq.gz, name1_L002_XYZ_R2_001.fastq.gz")
    print("- name1_L002_XYZ_R1_002.fastq.gz, name1_L002_XYZ_R2_002.fastq.gz\n")
    
    print("\t 1.By adding the option --include_all in all modules, XICRA will identify four samples:")
   
    print("\t - Sample 1: name1_L001_XYZ_R1_001, name1_L001_XYZ_R2_001")
    print("\t - Sample 2: name1_L001_XYZ_R1_002, name1_L001_XYZ_R2_002")
    print("\t - Sample 3: name1_L002_XYZ_R1_001, name1_L002_XYZ_R2_001")
    print("\t - Sample 4: name1_L002_XYZ_R1_002, name1_L002_XYZ_R2_002")
    print (colored('\t ** Remember to use option --include_all within each module. **\n', 'yellow'))
    
    print("\t 1.By adding the options --include_all --merge_Reads within module prep, XICRA will identify two samples:")
    print("\t - Sample 1: name1_L001_XYZ_R1, name1_L001_XYZ_R2")
    print("\t - Sample 2: name1_L002_XYZ_R1, name1_L002_XYZ_R2")
    print (colored('\t ** Remember to use option --include_all within each module. **\n', 'yellow'))
    
def project_help():
    return ()

###############
def multiqc_help():
    return ()

###############
def print_help_adapters():
    print("Check the XICRA documentation of the module: ", "\033[4m"+ "https://xicra.readthedocs.io/en/latest/user_guide/modules/trimm.html"+"\033[0m")
    return()

###############
def help_join_reads():
    return ()

###############
def help_miRNA():
    print("Check the XICRA documentation of the module: ", "\033[4m"+ "https://xicra.readthedocs.io/en/latest/user_guide/modules/miRNA.html"+"\033[0m")
    return()

###############
def help_tRNA():
    return ()

###############
def help_piRNA():
    return ()
