import os,sys,re
from python_bymyself import read_fasta


NBKK_pro=sys.argv[1]
Niben_pro=sys.argv[2]
blastp_NBKK2Niben=sys.argv[3]
blastp_Niben2NBKK=sys.argv[4]
tblastn_NBKK2Niben=sys.argv[5] ##tblastn to Niben genome
tblastn_Niben2NBKK=sys.argv[6]
NBKK_gff=sys.argv[7]
Niben_gff=sys.argv[8]
NBKK_DRAGO=sys.argv[9]
Niben_DRAGO=sys.argv[10]
Niben_all=sys.argv[11]

def dragodic(file):
    dic_out={}
    with open(file,'r') as file_in:
        for line in file_in:
            geneid=line.split()[0]
            if re.match(r'Nbe',geneid):
                geneid=re.sub(r'\.\d+$','',geneid)
            dic_out[geneid]=line.split()[1]
    return dic_out

def retrive_blastp(file,fastalen):
    dic_out={}
    exists={}
    with open(file,'r') as file_in:
        for line in file_in:
            queryid=line.split()[0]
            sbjectid=line.split()[1]
            if queryid in exists.keys():
                pass
            else:
                matchratio=line.split()[2]
                matchlen=line.split()[3]
                querymatch=str(abs(int(line.split()[6])-int(line.split()[7]))+1)
                sbjectmatch=str(abs(int(line.split()[8])-int(line.split()[9]))+1)
                dic_out[queryid]=sbjectid+'\t'+matchratio+'\t'+querymatch+'\t'+sbjectmatch
            
            exists[queryid]=1
    return dic_out
def retrive_tblastn(file):
    dic_out={}
    exists={}
    with open(file,'r') as file_in:
        for line in file_in:
            queryid=line.split()[0]
            sbjectid=line.split()[1]
            if queryid in exists.keys():
                pass
            else:
                matchratio=float(line.split()[2])
                if int(line.split()[8])<int(line.split()[9]):
                    dic_out[queryid]=sbjectid+':'+line.split()[8]+'-'+line.split()[9]+'+'
                else:
                    dic_out[queryid]=sbjectid+':'+line.split()[9]+'-'+line.split()[8]+'-'
            
            exists[queryid]=1
    return dic_out
def retrive_gff(file):
    dic_out={}
    with open(file,'r') as file_in:
        for line in file_in:
            if re.search(r'\tgene\t',line):
                geneId=re.search(r'ID=(\S+)',line).group(1)
                geneId=re.sub(r'\;','',geneId)
                dic_out[geneId]=line.split()[0]+':'+line.split()[3]+'-'+line.split()[4]+line.split()[6]
    return dic_out
dic_NBKKfalen=read_fasta.fasta_len(NBKK_pro)
dic_Nibenfalen=read_fasta.fasta_len(Niben_all)
NBKKfaid=read_fasta.fasta_id_list(NBKK_pro)
Nibenfaid=read_fasta.fasta_id_list(Niben_pro)
dic_blastpKK2Niben=retrive_blastp(blastp_NBKK2Niben,dic_NBKKfalen)
dic_blastpNiben2KK=retrive_blastp(blastp_Niben2NBKK,dic_Nibenfalen)
dic_NBKKgene2loci=retrive_gff(NBKK_gff)
dic_Nibengene2loci=retrive_gff(Niben_gff)
dic_NBKK2Nibengenome=retrive_tblastn(tblastn_NBKK2Niben)
dic_Niben2NBKKgenome=retrive_tblastn(tblastn_Niben2NBKK)
NBKKgene2name=dragodic(NBKK_DRAGO)
Nibengene2name=dragodic(Niben_DRAGO)
#Gene Type	NBKK Gene ID	NBHZ Gene ID	Coverage	Percent Identity	NBKK genomic location	NBHZ genomic location
print('NBKK Gene ID\tNBKK Gene type\tNBHZ Gene ID\tNBHZ Gene type\tNBKK coverage\tNBHZ coverage\tPercent Identity\tNBKK genomic location\tNBHZ genomic location')
removeNibenid={}
dic_final={}
for eachid in sorted(NBKKfaid):
    if eachid in dic_blastpKK2Niben.keys():
        matchNibenid=dic_blastpKK2Niben[eachid].split('\t')[0]
        KKgeneid=re.sub(r'\.\d+$','',eachid)
        Nibengeneid=re.sub(r'\.\d+$','',matchNibenid)
        kkgeneout=eachid
        if eachid in NBKKgene2name.keys():
            kkgeneout=kkgeneout+'\t'+NBKKgene2name[eachid]
        else:
            kkgeneout=kkgeneout+'\t-'
        nbgeneout=matchNibenid
        if Nibengeneid in Nibengene2name.keys():
            nbgeneout=nbgeneout+'\t'+Nibengene2name[Nibengeneid]
        else:
            nbgeneout=nbgeneout+'\t-'
        #if matchNibenid in dic_blastpNiben2KK.keys():
            #reverseNBKKid=dic_blastpNiben2KK[matchNibenid].split('\t')[0]
        removeNibenid[matchNibenid]='Y'
        matchratio=dic_blastpKK2Niben[eachid].split('\t')[1]
        NBKKmatchlen=dic_blastpKK2Niben[eachid].split('\t')[2]
        Nibenmatchlen=dic_blastpKK2Niben[eachid].split('\t')[3]
        NBKKcover=round(float(NBKKmatchlen)/dic_NBKKfalen[eachid]*100,2)
        Nibencover=round(float(Nibenmatchlen)/dic_Nibenfalen[matchNibenid]*100,2)
        dic_final.setdefault(KKgeneid,[])
        dic_final[KKgeneid].append(kkgeneout+'\t'+nbgeneout+'\t'+str(NBKKcover)+'\t'+str(Nibencover)+'\t'+matchratio+'\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_Nibengene2loci[Nibengeneid])
        #print(kkgeneout+'\t'+nbgeneout+'\t'+str(NBKKcover)+'\t'+str(Nibencover)+'\t'+matchratio+'\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_Nibengene2loci[Nibengeneid])
        
for eachone in Nibenfaid: 
    if eachone in removeNibenid.keys():
        pass
    else:
        matchKKid=dic_blastpNiben2KK[eachone].split('\t')[0]
        KKgeneid=re.sub(r'\.\d+$','',matchKKid)
        Nibengeneid=re.sub(r'\.\d+$','',eachone)
        kkgeneout=matchKKid
        if matchKKid in NBKKgene2name.keys():
            kkgeneout=kkgeneout+'\t'+NBKKgene2name[matchKKid]
        else:
            kkgeneout=kkgeneout+'\t-'
        nbgeneout=eachone
        if Nibengeneid in Nibengene2name.keys():
            nbgeneout=nbgeneout+'\t'+Nibengene2name[Nibengeneid]
        else:
            nbgeneout=nbgeneout+'\t-'
        #if matchNibenid in dic_blastpNiben2KK.keys():
            #reverseNBKKid=dic_blastpNiben2KK[matchNibenid].split('\t')[0]
        matchratio=dic_blastpNiben2KK[eachone].split('\t')[1]
        NBKKmatchlen=dic_blastpNiben2KK[eachone].split('\t')[3]
        Nibenmatchlen=dic_blastpNiben2KK[eachone].split('\t')[2]
        NBKKcover=round(float(NBKKmatchlen)/dic_NBKKfalen[KKgeneid]*100,2)
        Nibencover=round(float(Nibenmatchlen)/dic_Nibenfalen[eachone]*100,2)
        dic_final.setdefault(KKgeneid,[])
        dic_final[KKgeneid].append(kkgeneout+'\t'+nbgeneout+'\t'+str(NBKKcover)+'\t'+str(Nibencover)+'\t'+matchratio+'\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_Nibengene2loci[Nibengeneid])
        #print(kkgeneout+'\t'+nbgeneout+'\t'+str(NBKKcover)+'\t'+str(Nibencover)+'\t'+matchratio+'\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_Nibengene2loci[Nibengeneid])
existsid={}                
for eachid in sorted(list(dic_final.keys())):
    for eachp in dic_final[eachid]:
        thiskey=eachp.split()[0]+'\t'+re.sub(r'\.\d$','',eachp.split()[2])
        if thiskey in existsid.keys():
            pass
        else:
            print(eachp)
        existsid[thiskey]='Y'
        
        
'''        
        KKgeneid=re.sub(r'\.\d+$','',eachid)
        Nibengeneid=re.sub(r'\.\d+$','',matchNibenid)
        kkgeneout=eachid
        if eachid in NBKKgene2name.keys():
            kkgeneout=kkgeneout+':'+NBKKgene2name[eachid]
        nbgeneout=matchNibenid
        if matchNibenid in Nibengene2name.keys():
            nbgeneout=nbgeneout+':'+Nibengene2name[matchNibenid]
        if matchNibenid in dic_blastpNiben2KK.keys():
            if dic_blastpNiben2KK[matchNibenid]==eachid:
                print('P2P'+'\t'+kkgeneout+'\t'+nbgeneout+'\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_Nibengene2loci[Nibengeneid])
                Nibenfaid.remove(matchNibenid)
            else:
                print('NBKK'+'\t'+kkgeneout+'\t'+nbgeneout+'\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_Nibengene2loci[Nibengeneid])
        else:
            print('NBKK'+'\t'+kkgeneout+'\t'+nbgeneout+'\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_Nibengene2loci[Nibengeneid])
    elif eachid in dic_NBKK2Nibengenome.keys():
        KKgeneid=re.sub(r'\.\d+$','',eachid)
        kkgeneout=eachid
        if eachid in NBKKgene2name.keys():
            kkgeneout=kkgeneout+':'+NBKKgene2name[eachid]
        print('NBKK\t'+kkgeneout+'\tno homolog\t'+dic_NBKKgene2loci[KKgeneid]+'\t'+dic_NBKK2Nibengenome[eachid])
    else:
        kkgeneout=eachid
        if eachid in NBKKgene2name.keys():
            kkgeneout=kkgeneout+':'+NBKKgene2name[eachid]
        print('O2K\t'+kkgeneout+'\tonly in NBKK')
for eachid in sorted(Nibenfaid):
    nbgeneout=eachid
    if eachid in Nibengene2name.keys():
        nbgeneout=nbgeneout+':'+Nibengene2name[eachid] 
    if eachid in dic_Niben2NBKKgenome.keys():
        print('Niben\tno homology\t'+nbgeneout+'\t'+dic_Niben2NBKKgenome[eachid]+'\t'+dic_Nibengene2loci[re.sub(r'\.\d+$','',eachid)])
    else:
        print('O2N\t'+nbgeneout)
'''        



