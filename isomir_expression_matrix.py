#usr/bin/env python

"""
    Author: Antonio Luna de Haro
    Internship in Bioinformatics (Lauro Sumoy)
    
    input --> un prefix per a la sample que vols. ha de tenir fet els isomirs.gtf 
    output --> escriu taulas de expressio individuals per a cada sample 
"""

from sys import argv
import os


def getfiles_isomir(samples_prefix):
    path_s = '/imppc/labs/lslab/aluna/isomiRs/single/' #'BST_10_L001/Stats'
    path_p = '/imppc/labs/lslab/aluna/isomiRs/paired/' #'BST_10_L001/Stats'
    files_s = os.listdir(path_s)
    sample_s_list = []
    for singles in files_s:
        if singles.startswith(samples_prefix): 
           sample_s_list.append(path_s + singles + '/' + singles + '.gtf' )
    
    files_p = os.listdir(path_p)
    sample_p_list = []
    for paireds in files_p:
        if paireds.startswith(samples_prefix): 
           sample_p_list.append(path_p + paireds + '/' + paireds + '.gtf')
            
    isomirs_s = list(set(sample_s_list))
    isomirs_p = list(set(sample_p_list))

    return sorted(isomirs_s), sorted(isomirs_p)
    #return isomir_dict    

def parse_gtf(gtffile):
    gtfile = open(gtffile)
    text = gtfile.read()
    lines = text.splitlines()
    sample_dict = {}
    for line in lines:
        if not line.startswith('#'):
            name = line.split('\t')[0]
            name = line.split('\t')[-1].split(';')[2].split()[-1]
            variant = line.split('\t')[-1].split(';')[4].split()[-1]
            variant = '-' + variant
            if variant == '-NA':
                variant = ''
            expression = int(line.split('\t')[-1].split(';')[6].split()[-1])
            namevariant = name + variant
            if namevariant in sample_dict.keys():
                old = sample_dict[namevariant]
                new = old + expression
                sample_dict[namevariant] = new
            else:
                sample_dict[namevariant] = expression
            #sample_dict[name+variant]= {}
    return sample_dict
            
            
  
    
    
def make_table(dictionary, filename):
    fil = open(filename, 'w')    
    isomirs = dictionary.keys()
    fil.write('Isomir\t')
    fil.write(filename.partition('_expression')[0])
    fil.write('\n')
    for sample in isomirs:
        fil.write(sample)
        fil.write('\t')
        fil.write(str(dictionary[sample]))
        fil.write('\n')
    fil.close()    
    

if __name__ == "__main__":
    #sample_dict = parse_gtf('test.gtf')
    #print sample_dict

#    testfile = argv[1]
#    dicti = parse_gtf(testfile)
#    make_table(dicti, 'test_ex_ma.txt')
    prefix1 = argv[1]
    
    #prefix2 = argv[2]
    singles_list, paired_list = getfiles_isomir(prefix1)
    #singles_list2, paired_list2 = getfiles_isomir(prefix2)
    filenames = []
    print 'Starting with sinles'
    for gtffile in singles_list:
        sample = gtffile.rpartition('/')[-1][:-4]
        sample_dict = parse_gtf(gtffile)
        filename = '/imppc/labs/lslab/aluna/LastTry/mirnaisomir/single/' + sample + '_expression_matrix.txt'
        filenames.append(filename)
        make_table(sample_dict, filename)
        print sample, ' expression matrix created'
    
    print 'Singles, done, now paired'
    for gtffile in paired_list:
        sample = gtffile.rpartition('/')[-1][:-4]
        sample_dict = parse_gtf(gtffile)
        filename = '/imppc/labs/lslab/aluna/LastTry/mirnaisomir/paired/' + sample + '_expression_matrix.txt'
        filenames.append(filename)
        make_table(sample_dict, filename)
        print sample, ' expression matrix created'
        
    
    #print ",".join(filenames)

    
    
    