###############   
def parse_gtf(gtffile, filename, sample, tag):

    gtfile = open(gtffile)
    text = gtfile.read()
    lines = text.splitlines()
    sample_dict = {}
    
    ## Open file
    fil = open(filename, 'w')    
    string2write = 'type\tsample_name\tparent\tname\tvariant\tUID\tseq\texpression\n'
    fil.write(string2write)
    for line in lines:
        if not line.startswith('#'):        
            #print ('## ',line)
            seq = line.split('\t')[-1].split(';')[0].split("Read=")[-1]
            ident = line.split('\t')[0]
            name = line.split('\t')[-1].split(';')[2].split("Name=")[-1]
            UID = line.split('\t')[-1].split(';')[1].split("UID=")[-1]
            parent = line.split('\t')[-1].split(';')[3].split("Parent=")[-1]
            variant = line.split('\t')[-1].split(';')[4].split("Variant=")[-1]
            expression = str(line.split('\t')[-1].split(';')[6].split("Expression=")[-1])
            string2write = tag + '\t' + sample + '\t' + ident + '\t' + name + '\t' + variant + '\t' + UID + '\t' + seq + '\t' + expression + '\n'
            fil.write(string2write)
    
    fil.close()      
###############