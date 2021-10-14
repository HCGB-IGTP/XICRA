

























#===============================================================================
# def MINTmap(reads, folder, file_name, num_threads):
#     MINTmap = config['EXECUTABLES']['MINTmap_folder'] + 'MINTmap.pl'    
#     MINTmap_table = config['EXECUTABLES']['MINTmap_folder'] + 'LookupTable.tRFs.MINTmap_v1.txt'
#     MINTmap_tRNAseq = config['EXECUTABLES']['MINTmap_folder'] + 'tRNAspace.Spliced.Sequences.MINTmap_v1.fa'
#     MINTmap_tRF = config['EXECUTABLES']['MINTmap_folder'] + 'OtherAnnotations.MINTmap_v1.txt'    
#     MINTmap_MINTplates = config['EXECUTABLES']['MINTmap_folder'] + 'MINTplates/'
#     results = []
#     command2sent = []
#     
#     ## open file
#     output_file = open(file_name, 'a')
#     output_file.write("\nMINTmap:\n")
#     
#     for jread in reads:    
#         for prefix in prefix_list:        
#             if paired_end:
#                 sample_search = re.search(r"(%s)\_(\d{1,2})\_(.*)" % prefix, jread)
#             else:
#                 name_sample = os.path.basename(jread)
#                 name_dir = os.path.dirname(jread)
#                 sample_search = re.search(r"(.*)\_trimmed\.fastq", name_sample)
#                 
#             if sample_search:
#                 if paired_end:
#                     outdir = sample_search.group(1) + "_" + sample_search.group(2)
#                 else:
#                     outdir = sample_search.group(1)
# 
#                 sample_folder =  folder + '/' + outdir + '/'
#                 results.append(sample_folder)
#                 logfile = sample_folder + outdir + '_logfile.txt'
#                 if (os.path.isdir(sample_folder)):
#                     print ('\tMINTmap analysis for sample %s already exists' %outdir)         
#                 else:
#                     #MINTmap.pl -f trimmedfastqfile [-p outputprefix] [-l lookuptable] [-s tRNAsequences] [-o tRFtypes] [-d customRPM] [-a assembly] [-j MINTplatesPath] [-h]
#                     fol = functions.create_subfolder(outdir, folder)
#                     cmd = 'perl '+ MINTmap + ' -f %s -p %s -l %s -s %s -o %s -j %s > %s' %(jread, sample_folder + outdir, MINTmap_table, MINTmap_tRNAseq, MINTmap_tRF, MINTmap_MINTplates, logfile) 
#                     # get command
#                     command2sent.append(cmd)
#                     # print into file
#                     output_file.write(cmd)
#                     output_file.write('\n')
# 
#     #sent commands on threads            
#     output_file.close()
#     
#     #print ("Commands:")
#     #print (len(command2sent))
#     functions.sender(command2sent, num_threads)
# 
#     results = set(results) ## BUG: if single-end option, it sends as many as prefixes each command 
#     #print ("results:")
#     #print (len(results))
#     
#     return results        
# ###############
# 
# ###############
# def tRFs_analysis(path, count, reads, time_partial, output_file, num_threads):
#     ##############################
#     ####### Step: sRNAbench ######
#     ##############################
#     print ("\n+ Run MINTmap: ")
#     name_MINTmap_folder = str(count) + '.1.tRFs_MINTmap'
#     MINTmap_folder = functions.create_subfolder(name_MINTmap_folder, path)
#     results = MINTmap(reads, MINTmap_folder, output_file, num_threads)
#     
#     print ("\n+ Get MINTmap matrix: ")
#     name_MINTmap_matrix = str(count) + '.2.tRFs_matrix'
#     MINTmap_matrix_folder = functions.create_subfolder(name_MINTmap_matrix, path)
# 
#     for folder in results:
#         files = os.listdir(folder)
#         for item in files:
#             if 'countsmeta' in item:
#                 continue
#             if item.endswith('html'):
#                 continue
#             if 'ambigu' in item:
#                 parse_tRF(folder, item, MINTmap_matrix_folder, 'ambiguous')        
#             elif 'exclu' in item:
#                 parse_tRF(folder, item, MINTmap_matrix_folder, 'exclusive')        
#     
#     ## functions.timestamp
#     time_partial = functions.timestamp(time_partial)
#     
# ###############
# def parse_tRF(path, fileGiven, matrix_folder, ident):
#     pathFile = path + '/' + fileGiven
#     sample_search = re.search('(.*)\-MINTmap_v1.*', fileGiven)    
#     if sample_search:
#         sample_name = sample_search.group(1)
#         #sample_folder =  matrix_folder + '/' + sample_name
#         skip = 0
#         tsv_file = matrix_folder + '/' + sample_name + '_' + ident + '.tsv'
#         if os.path.isfile(tsv_file):
#             print ('\tMatrix for ', sample_name , ' (' + ident + ') is already generated')
#             skip = 1
#         if skip == 0:
#              ## Open file
#             fil = open(tsv_file, 'w')
#             string2write = 'type\tsample_name\tident\tname\tvariant\tUID\tseq\texpression\n'
#             fil.write(string2write)
#             ## Read file
#             expression_file = open(pathFile)
#             expression_text = expression_file.read()
#             expression_lines = expression_text.splitlines()
#             for line in expression_lines:
#                 if not line.startswith('MINTbase'):
#                     UID = line.split('\t')[0]
#                     seq = line.split('\t')[1]
#                     variant = line.split('\t')[2]
#                     expression = line.split('\t')[3]
#                     tRNA_name = line.split('\t')[-1].split(',')[0]
# 
#                     # tRF-31-87R8WP9N1EWJ0    TCCCTGGTGGTCTAGTGGTTAGGATTCGGCG    5'-tRF    921    7026.67    452.60    na    trna77_GluCTC_6_+_28949976_28950047@1.31.31, trna80_GluCTC_1_-_161417018_161417089@1.31.31
#                     tRNA_search = re.search(r"trna.{1,3}\_(.{6})\_(.{1,2})\_.*", tRNA_name)
#                     tRNA_family = 'na'
#                     if tRNA_search:
#                         tRNA_family = tRNA_search.group(1)
#                         if (tRNA_search.group(2) == 'MT'):
#                             tRNA_family = tRNA_family + '_MT'
#                         
#                     
#                     string2write = 'tRFs\t'+ sample_name + '\t' + ident + '\t' + tRNA_family +'\t' + variant +'\t' + UID + '\t' + seq + '\t' + expression + '\n'
#                     fil.write(string2write)
# 
#             fil.close()
# ###############
#===============================================================================
