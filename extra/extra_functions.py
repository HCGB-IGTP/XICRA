def select_other_samples (project, list_samples, samples_prefix, mode, extensions, exclude=False, Debug=False):

    ## init dataframe
    name_columns = ("sample", "dirname", "name", "ext", "tag")

    ## initiate dataframe
    df_samples = pd.DataFrame(columns=name_columns)

    #Get all files in the folder "path_to_samples"    
    sample_list = []
    for names in samples_prefix:
        for path_file in list_samples:    
            f = os.path.basename(path_file)
            dirN = os.path.dirname(path_file)
            #samplename_search = re.search(r"(%s).*" % names, f)
            samplename_search = re.search(r"(%s).*" % names, path_file)
            
            enter = ""
            if samplename_search:
                if (exclude): ## exclude==True
                    enter = False
                else: ## exclude==True
                    enter = True
            else:
                if (exclude): ## exclude==True
                    enter = True
                else: ## exclude==True
                    enter = False
                    
            if enter:
                
                ## project mode:
                if project:
                    if mode == 'annot':
                        #### /path/to/folder/annot/name.faa
                        for ext in extensions:
                            f_search = re.search(r".*\/%s\/(.*)\.%s$" %(mode, ext), path_file)
                            if f_search:
                                file_name = f_search.group(1) 
                                df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, ext, mode]    

                    elif mode== 'assembly':
                        #### name_assembly.faa
                        f_search = re.search(r"(.*)\_%s\.%s$" %(mode, extensions), f)
                        if f_search:
                            file_name = f_search.group(1) 
                            df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, extensions, mode]    

                    elif mode== 'mash':
                        #### name.sig
                        f_search = re.search(r".*\/%s\/(.*)\.%s$" %(mode, extensions[0]), path_file)
                        if f_search:
                            file_name = f_search.group(1) 
                            df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, extensions[0], mode]    

                    else:
                        f_search = re.search(r".*\/(.*)\/%s\/(.*)\_summary\.%s$" %(mode, extensions[0]), path_file)
                        if f_search:
                            ### get information
                            if mode == 'profile':
                                name = f_search.group(1)
                                db_name = f_search.group(2).split('_')[-1]
                                if not name.startswith('report'):
                                    df_samples.loc[len(df_samples)] = [path_file, dirN, name, db_name, mode]    

                            elif mode == 'ident':
                                name = f_search.group(1)
                                df_samples.loc[len(df_samples)] = [path_file, dirN, name, 'csv', mode]    

                ## detached mode
                else:
                    if f.endswith(extensions):
                        file_name, ext = os.path.splitext(f)
                        df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, db_name, mode]    
                    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: df_samples **", 'yellow'))
        print (df_samples)
    
    ##
    number_samples = df_samples.index.size
    if (number_samples == 0):
        print (colored("\n**ERROR: No samples were retrieved for this option. Continue processing...\n",'red'))
        return (df_samples)
    print (colored("\t" + str(number_samples) + " samples selected from the input provided...", 'yellow'))

    return (df_samples)
