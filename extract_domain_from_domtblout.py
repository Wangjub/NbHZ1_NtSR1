import re,os,sys
from python_bymyself import read_fasta
domtblout=sys.argv[1]
protein_db=sys.argv[2]
prefix=sys.argv[3]
output=sys.argv[4]

dic_all={}
for eachfile in protein_db.split(','):
    each_fasta=read_fasta.read_fasta(eachfile)
    dic_all.update(each_fasta)

dic_domain={}    
outwrite=open(output,'w')
with open(domtblout,'r') as file_in:
    for line in file_in:
        if re.match(prefix,line):
            geneId=line.split()[3]
            dic_domain.setdefault(geneId,[])
            dic_domain[geneId].append(line.strip())
            
for eachid,eachunit in dic_domain.items():
    unitsort=sorted(eachunit,key=lambda x:-abs(int(x.split()[17])-int(x.split()[18])))
    start=int(unitsort[0].split()[17])-1
    end=int(unitsort[0].split()[18])-1
    thisseq=dic_all[eachid][start:end]
    outwrite.write('>'+eachid+'\n'+thisseq+'\n')
outwrite.close()

