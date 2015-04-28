#!/usr/bin/env python
# -*- coding: utf8 -*-

from Bio import SeqIO
from multiprocessing import Process, Manager
import os, sys, random, re


input_file = sys.argv[1]
list1=[]
spec_codes=[]

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
    list1.append(fasta.id)


output_fna=sys.argv[2]
output_bed=sys.argv[3]
type=sys.argv[4]
flank=sys.argv[5]
begin=sys.argv[6]
end=sys.argv[7]

for i in list1:
	toto = i.split('_')
	spec_codes.append(toto[1])

myset = set(spec_codes)
print myset

	
path_file=open('/usr/local/bioinfo/galaxy_dev/galaxy_dist/tool-data/gff2fasta.loc','r')
spec=path_file.readlines()
path_file.close()
bases_gff=[]
bases_fna=[]
for line in spec :
	infos=line.split('\t')
	if infos[0] in spec_codes:
		bases_gff.append(infos[2].rstrip())
		bases_fna.append(infos[3].rstrip())
		print(infos[2].rstrip()+"\t"+infos[3].rstrip())

curdir=os.getcwd()
    	
numb=str(random.randint(1,100000))
tmpfoldname=curdir+"/"+numb
os.system("mkdir "+tmpfoldname)
os.system("mkdir "+tmpfoldname+"/logs")
os.system("mkdir "+tmpfoldname+"/outputs")


i=1
for path in bases_gff:
	os.system("ln -s "+ path +" "+tmpfoldname+"/gff"+str(i))
	#print path
	i=i+1
	
j=1	
for path in bases_fna:
	os.system("ln -s "+ path +" "+tmpfoldname+"/fna"+str(i))
	#print path
	j=j+1


jobarray="#!/bin/bash\n#$ -N gff2fnaclust\n#$ -wd "+tmpfoldname+"/\n#$ -e "+tmpfoldname+"/logs/ \n#$ -o "+tmpfoldname+"/logs/\n#$ -q bioinfo.q\n#$ -t 1-"+str(i-1)+"\n#$ -tc "+str(i-1)+"\n#$ -S /bin/bash \n#$ -b y\n#$ -V\nperl /usr/local/bioinfo/galaxy/galaxy_dist/tools/gff2fna/multi_gff2fna.pl ./gff${SGE_TASK_ID} ./fna${SGE_TASK_ID} "+output_fna+" "+output_bed+" "+input_file+" "+type+" "+flank+" "+begin+" "+end
array_file=open(tmpfoldname+'/jobs.sge','w')
array_file.write(jobarray)
array_file.close()


os.system("qsub -sync y "+tmpfoldname+"/jobs.sge")
os.chdir(tmpfoldname)

os.chdir("../")
os.system("rm -rf "+tmpfoldname)