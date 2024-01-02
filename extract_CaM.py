from Bio.Blast import NCBIXML
import re,sys
from python_bymyself import read_fasta
blastResult=sys.argv[1]
protein_db=sys.argv[2]
all_blast=sys.argv[3]
query_fa=sys.argv[4]
dic_fasta=read_fasta.read_fasta(protein_db)

dic_SR2other={}
exists_id=[]
number=0
with open(all_blast,'r') as file_in:
    for line in file_in:
        geneid=line.split()[0]
        matchid=line.split()[1]
        matchspecies=''
        number+=1
        #print(number,end='\r')
        #print(line,end='')
        if re.match(r'Soly',matchid):
            matchspecies='SLY'
        elif re.match(r'LOC',matchid):
            matchspecies='OSA'
        elif re.match(r'AT',matchid):
            matchspecies='ATH'
        
        if geneid+'\t'+matchspecies in exists_id:
            pass
        else:
            matchid2=re.sub(r'\.\S+$','',matchid)
            dic_SR2other.setdefault(geneid,[])
            dic_SR2other[geneid].append(matchid2)
        exists_id.append(geneid+'\t'+matchspecies)

#print(dic_SR2other)
queryList=[]
with open(query_fa,'r') as file_in:
    for line in file_in:
        if re.match(r'>',line):
            queryList.append(re.sub(r'\.\S+$','',re.match(r'>(\S+)',line).group(1)))


blastout=open(blastResult)
dic_gene={}
blast_records=NCBIXML.parse(blastout)
for blast_record in blast_records:
    for alignments in blast_record.alignments:
        for hsp in alignments.hsps:
            if hsp.align_length>blast_record.query_length*0.85 and hsp.identities/hsp.align_length>0.9 and alignments.length==149:
                #geneid=re.sub(r'\.t\d+','',alignments.hit_id)
                geneid=alignments.hit_id
                homologList=dic_SR2other[geneid]
                if not set(queryList).isdisjoint(set(homologList)):
                    dic_gene[geneid]=1
blastout.close()
for eachid in dic_gene.keys():
    seq2=re.sub(r'[\.\*]{1,}$','',dic_fasta[eachid])
    print('>'+eachid+'\n'+seq2)


