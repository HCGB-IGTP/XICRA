#usr/bin/env python

"""
    Author v0.1: Antonio Luna de Haro
    Author v1: Jose F Sanchez-Herrero

    This scripts runs sRNAtoolbox for a selection of RNAseq samples
!!!!    Make sure to activate miRTop! (activate-mirtop-0.3.11a0)    !!!!
"""
## useful imports
import sys
from sys import argv
import time
from datetime import datetime
import subprocess
import os
import re
import configparser


#####################
#### functions ######
#####################

def select_samples (samples_prefix = "", path_to_samples = ""):

    #Debug: print samples_prefix
    #Debug: print path_to_samples
    print 'Selecting samples from:', path_to_samples
    
    #Get all files in the folder "path_to_samples"    
    files = os.listdir(path_to_samples)
    sample_list = []
    for fastq in files:
        if fastq.startswith(samples_prefix) and 'merged' not in fastq: 
            if fastq.endswith('.gz'):
                sample_list.append(fastq[:-3]) 
            elif fastq.endswith('fastq'):
                sample_list.append(fastq)
            else:
                print "** ", fastq, 'is a file that is neither in fastq.gz or .fastq format, so it is not included'
    
    non_duplicate_samples = list(set(sample_list))
    number_samples = len(non_duplicate_samples)
    print '\nPipeline will run with all', samples_prefix, 'samples (%i samples)' %number_samples, '\n'
    return sorted(non_duplicate_samples)


def gunziping (non_duplicate_samples, reverse = False, path_to_samples = "/imppc/labs/lslab/share/isomiR_project/mirXPlore/FASTQ"):
    path_to_samples = "/imppc/labs/lslab/aluna/Ilumina_NGS/150120_VHIO_data_FIS_SERUM_miRNAs/fastqR"
    if reverse:
        path_to_samples = path_to_samples + '2'
    else:
        path_to_samples = path_to_samples + '1'
    
    final_sample_list = []     
    for sample in non_duplicate_samples:
        sample = path_to_samples + sample
        if not os.path.isfile(sample):
            sample = sample + '.gz'
            print 'Using gunzip to -->', sample
            subprocess.check_output(['gunzip','-k', sample])
            sample = sample[:-3]
        final_sample_list.append(sample)
    return sorted(final_sample_list)
    
    
    
def one_file_per_sample(final_sample_list):
    grouped_subsamples = []
    bigfile_list = []
    for samplex in final_sample_list:
    	print final_sample_list
        if samplex not in grouped_subsamples:
            subsamples = []
            samplename_search = re.search('(.*\/)([a-zA-Z]{2,3})(\d{5})(\d{2})(.*)', samplex)
            if samplename_search:
                path = samplename_search.group(1)
                sample = samplename_search.group(2)
                comonpart = path + sample
                print comonpart
                for sampley in final_sample_list:
                    if comonpart in sampley:
                        subsamples.append(sampley)
                        grouped_subsamples.append(sampley)
            bigfilepath = comonpart + '_merged.fastq'
            bigfilename = sample + '_merged.fastq'
            bigfile_list.append(bigfilepath)
            #print sorted(subsamples)
            if not os.path.isfile(bigfilepath) or os.stat(bigfilepath).st_size == 0:
                partsofsample = ' '.join(sorted(subsamples))
                cmd = 'cat %s > %s' %(partsofsample, bigfilepath)
                try:
                    print '%s created' %bigfilename
                    subprocess.check_output(cmd, shell = True)
                    
                    
                except subprocess.CalledProcessError as err:
                    print err.output
                    sys.exit()
            else:
                print 'Sample %s is already merged' %sample
    print 'There is %i samples after merging' %len(bigfile_list)
    return bigfile_list
    
   
def cutadapt (list_R1, list_R2):
    out_path = '/imppc/labs/lslab/aluna/Ilumina_NGS/150120_VHIO_data_FIS_SERUM_miRNAs/cutadapt_out/'
    trimed_R1 = []
    trimed_R2 = [] 
    for file_R1 in list_R1:
        
        sampleR1_search = re.search('(.*\/)([a-zA-Z]{2,3})(\d{5})(\d{2})(.*)', file_R1)
        if sampleR1_search:
            o_param = out_path + sampleR1_search.group(2) + '_' + sampleR1_search.group(4) + '_1.fastq'
            p_param = out_path + sampleR1_search.group(2) + '_' + sampleR1_search.group(4) + '_2.fastq'
        
            for file_R2 in list_R2:
                
                sampleR2_search = re.search('(.*\/)([a-zA-Z]{2,3})(\d{5})(\d{2})(.*)', file_R2)
                if sampleR2_search:
                    
                    if sampleR1_search.group(2) == sampleR2_search.group(2) and sampleR1_search.group(4) == sampleR2_search.group(4):
                        trimed_R1.append(o_param)
                        trimed_R2.append(p_param)
                        if not (os.path.isfile(o_param) or os.path.isfile(p_param)):
                            cmd = 'cutadapt -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC -o %s -p %s %s %s' %(o_param, p_param, file_R1, file_R2)
                            print 'The following cmd is being executed at the shell: \n', cmd
                            try:
                                print 'Trimming for %s sample' %(sampleR1_search.group(2) + '_' + sampleR1_search.group(4))
                                subprocess.check_output(cmd, shell = True)
                                print 'Adapters trimmed for the sample %s' %(sampleR1_search.group(2) + '_' + sampleR1_search.group(4))
                            except subprocess.CalledProcessError as err:
                                print err.output
                                sys.exit()                    
                        else:
                            print 'Sample %s is already trimmed' %(sampleR1_search.group(2) + '_' + sampleR1_search.group(4))
    return trimed_R1, trimed_R2

     
def fastqjoin (trimed_R1, trimed_R2, error_param = '0.1'):
    out_path = '/imppc/labs/lslab/aluna/Ilumina_NGS/150120_VHIO_data_FIS_SERUM_miRNAs/fastqjoin/'
    joined_reads = []
    for read1 in trimed_R1:
        sample1 = read1.rpartition('/')[2].split('.')[0][:-2]
        for read2 in trimed_R2:
            sample2 = read2.rpartition('/')[2].split('.')[0][:-2]
            if sample1 == sample2:
                joined_reads.append(out_path + sample1 + '_join.fastq')
                if not (os.path.isfile(out_path + sample1 + '_join.fastq')):
                    
                    cmd = '/soft/bio/ea-utils-expressionanalysis-git/bin/fastq-join -p %s %s %s -o %s' %(error_param, read1, read2, out_path + sample1 + '_%.fastq')
                    print 'The following cmd is being executed at the shell: \n', cmd
                    
                    
                    try:
                        print 'Sample %s is joining the paired-end reads' %sample1
                        subprocess.check_output(cmd, shell = True)
                        
                    except subprocess.CalledProcessError as err:
                        print err.output
                        sys.exit()
                else: 
                    print 'Sample %s is already joined' %sample1
    return joined_reads
    
    
    
def sRNAbench (joined_reads, paired=True):
    outpath = '/imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/out/'
    results = []
    
        
    for jread in joined_reads:
        if paired:
            outdir = jread.rpartition('/')[2].split('.')[0]
            finalpath = outpath + 'paired/' + outdir + '/'
        else:
            outdir = jread.rpartition('/')[2].split('.')[0][:-2]
            finalpath = outpath + 'single/' + outdir + '/'
        results.append(finalpath)
        print finalpath
        if not (os.path.isdir(finalpath)):
            cmd = 'java -jar /imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/exec/sRNAbench.jar dbPath=/imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/ input=%s output=%s microRNA=hsa  isoMiR=true plotLibs=true graphics=true plotMiR=true bedGraphMode=true writeGenomeDist=true chromosomeLevel=true chrMappingByLength=true ' %(jread, finalpath)
            print 'The following cmd is being executed at the shell: \n', cmd
            
            try:
                #print 'Searching isomiRs for %s' %outdir[:-5](.*\/)((([a-zA-Z]{2,3})_(C_)?(\d{1}|\d{2}))_[ACGT]{6}_(.{4})_R\d)_(.*)
                subprocess.check_output(cmd, shell = True)
                
            except subprocess.CalledProcessError as err:
                print err.output
                sys.exit()
        else: 
            print 'Sample %s has already been serached for isomiRs' %outdir[:-5]
    return results
    
def miRTop (results):
    outpath = '/imppc/labs/lslab/aluna/isomiRs/'
    gtfs = []
    for folder in results:
        outdir = '/'.join(folder.split('/')[-3:-1])
        print outdir
        gtfs.append(outpath + outdir + '/')
        if not (os.path.isdir(outpath + outdir)):
            cmd = 'mirtop gff --sps hsa --hairpin /imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/libs/hairpin.fa --gtf /imppc/labs/lslab/aluna/opt/hsa.gff3 --format srnabench -o %s  %s' %(outpath + outdir + '/', folder)

            print 'The following cmd is being executed at the shell: \n', cmd
            
            try:
                print 'Creating isomiRs gtf file for %s' %outdir
                subprocess.check_output(cmd, shell = True)
                
            except subprocess.CalledProcessError as err:
                print err.output
                sys.exit()
        else: 
            print 'Sample %s has already an isomiRs gtf file' %outdir
    return gtfs
    #return results
            
def miRTop_stats (gtfs):
    
    for gtf in gtfs:
        gtfname = gtf.split('/')[-2]
        gtfile = gtf + gtfname + '.gtf'
        outdir = gtf + 'Stats'
        if not (os.path.isdir(outdir)):
            cmd = 'mirtop stats -o %s %s' %(outdir, gtfile) 
            print 'The following cmd is being executed at the shell: \n', cmd
            
            try:
                print 'Creating isomiRs stats for %s' %gtfname
                subprocess.check_output(cmd, shell = True)
                
            except subprocess.CalledProcessError as err:
                print err.output
                sys.exit()
        else: 
            print 'Sample %s has already isomiRs stats' %gtfname

def getime (start_time):
    total_sec = time.time() - start_time
    m, s = divmod(int(total_sec), 60)
    h, m = divmod(m, 60)
    return h, m, s
    



if __name__ == "__main__":
    start_time_total = time.time()
  
    print "\n######## Starting Process ########"
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print date_time
    print "\n"
  
    ## get configuration
    configuration_path_file = os.path.abspath(argv[1])
    config = configparser.ConfigParser()
    config.read(configuration_path_file)

    print "------ Get configuration options ------"
    print "Reading configuration file: ", configuration_path_file
    print "\n"
    
    #Reading arguments from the comand line
    #if len(argv) > 2:
    #    
    #    s_prefix = argv[1]
    #    if argv[2] == 'single' or argv[2] == 'paired':
    #        type_reads = argv[2]
    #    else:
    #        print 'Second argument must be "single" or "paired". Using paired'
    #        type_reads = 'paired'
    #
    #else:
    #    print 'Using paired-end reads for BST samples (Two arguments must be provided if something else is wanted)'
    #    s_prefix = 'BST'
    #    type_reads = 'paired'
        
    ####### Step1: select_samples
    sample_listR1 = select_samples(samples_prefix="", path_to_samples=config["GENERAL"]["fastq_R1"])
    sample_listR2 = select_samples(samples_prefix="", path_to_samples=config["GENERAL"]["fastq_R2"])
    
    ####### Step2: merge_samples
    start_time_partial = time.time()
    mergedR1 = one_file_per_sample(sample_listR1)
    mergedR2 = one_file_per_sample(sample_listR2)
    h,m,s = getime(start_time_partial)
    print '\n---------------------------------------------------\n'
    print 'Each sample has been merged into a single file (Time spent: %i h %i min %i s)' %(int(h), int(m), int(s))
    print '\n---------------------------------------------------\n'

    exit()

    ####### Step3: cutadapt
    start_time_partial = time.time()
    trimed_R1, trimed_R2 = cutadapt(mergedR1, mergedR2)
    h,m,s = getime(start_time_partial)
    print '\n---------------------------------------------------\n'
    print 'Removed adapters for all samples (Time spent: %i h %i min %i s)' %(int(h), int(m), int(s))
    print '\n---------------------------------------------------\n'
    start_time_partial = time.time()
    joined_reads =fastqjoin(trimed_R1, trimed_R2)
    h,m,s = getime(start_time_partial)
    print '\n---------------------------------------------------\n'
    print 'Paired end reads converted into one (Time spent: %i h %i min %i s)' %(int(h), int(m), int(s))
    print '\n---------------------------------------------------\n'
    start_time_partial = time.time()
    if type_reads == 'paired':    
        results = sRNAbench(joined_reads)
    else:
        results = sRNAbench(trimed_R1, paired=False)
    h,m,s = getime(start_time_partial)
    print '\n---------------------------------------------------\n'
    print 'IsomiRs for each sample have been reported (Time spent: %i h %i min %i s)' %(int(h), int(m), int(s))
    print '\n---------------------------------------------------\n'
    start_time_partial = time.time()
    gtfs =miRTop(results)
    h,m,s = getime(start_time_partial)
    print '\n---------------------------------------------------\n'
    print 'Every sample has its own gtf file containing all isomiRs (Time spent: %i h %i min %i s)' %(int(h), int(m), int(s))
    print '\n---------------------------------------------------\n'
    start_time_partial = time.time()
    miRTop_stats(gtfs)
    h,m,s = getime(start_time_partial)
    ht,mt,st = getime(start_time_total)
    print '\n---------------------------------------------------\n'
    print 'EVERYTHING RAN WITHOUT PROBLEMS! EUREKA (Time spent: %i h %i min %i s)' %(int(h), int(m), int(s))
    print '(Total time spent: %i h %i min %i s)' %(int(ht), int(mt), int(st))
    print '\n---------------------------------------------------\n'
    
