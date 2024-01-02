import re,os,sys

from python_bymyself import read_fasta

input1=sys.argv[1]

dic_gene2trans={}
dic_id2seq=read_fasta.read_fasta(input1)
dic_id2len=read_fasta.fasta_len(input1)
with open(input1,'r') as file_in:
    for line in file_in:
        if re.match(r'>',line):
            transid=re.search(r'>(\S+)',line).group(1)
            geneid=re.sub(r'\.\S+$','',transid)
            dic_gene2trans.setdefault(geneid,[])
            dic_gene2trans[geneid].append(transid)


for eachgene in sorted(list(dic_gene2trans.keys())):
    transList=sorted(dic_gene2trans[eachgene],key=lambda x:-dic_id2len[x])
    transid2=transList[0]
    transseq2=dic_id2seq[transid2]
    transseq2=re.sub(r'[\.\*]{1,}$','',transseq2)
    print('>'+transList[0]+'\n'+transseq2)
