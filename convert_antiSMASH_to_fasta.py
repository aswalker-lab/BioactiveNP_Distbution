# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 14:01:09 2019

@author: Allison Walker
"""

from Bio import SeqIO
from os import listdir
from os.path import isfile, join, exists
#put parent directory for all antismash files ot be converted here
directory_name = "./"
def getAllFilesInDirectory(path, files):
    for f in listdir(path):
        if isfile(join(path, f)):
            #print join(path, f)
            if (".gbff" in f or ".gb" in f or ".gbk" in f) and ("region" in f or "cluster" in f) and "antismash5" in path:
                #check if output exists
                fullpath = join(path,f)
                if not exists(fullpath[0:fullpath.rfind(".")] + ".fasta"): 
                    #print(join(path, f))
                    files.append(join(path, f))
            #print files
        else:
            if "genometools" in f or "glimmer" in f:
                continue
            files = getAllFilesInDirectory(join(path, f), files)
    return files

files = []
files_to_convert = getAllFilesInDirectory(directory_name, files)
#print(files)

for f in files:
    record = SeqIO.read(open(f, 'rU'),"genbank")
    try:
        fasta_out = open(f[0:f.rfind(".")]+".fasta",'w')
        fasta_out.write(">")
        fasta_out.write(record.id + "|"+ record.description + "\n")
        for s in record.seq:
            fasta_out.write(s)
        fasta_out.write("\n")
        fasta_out.close()
    except:
        print("failed on: " + f)