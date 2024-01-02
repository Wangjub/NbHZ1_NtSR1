import re,os,sys
from python_bymyself import read_fasta

Blastf=sys.argv[1]
Nta_db=sys.argv[2]
Nbe_db=sys.argv[3]
Atpep=sys.argv[4]
Slpep=sys.argv[5]
Nr_Nt=sys.argv[6]
Nr_Nb=sys.argv[7]
Nt_gff=sys.argv[8]
Nb_gff=sys.argv[9]
out1=sys.argv[10]
out2=sys.argv[11]

def annotfile(filename):
    dic_out={}
    with open(filename,'r') as file_in:
        for line in file_in:
            if re.match(r'\S+',line):
                dic_out[line.split('\t')[0]]=line.split('\t')[1]
    return dic_out

def gff2loci(filename):
    dic_out={}
    with open(filename,'r') as file_in:
        for line in file_in:
            if re.search(r'\tmRNA\t',line):
                transid=re.search(r'ID=(\S+?)\;',line).group(1)
                dic_out[transid]=line.split('\t')[0]+':'+line.split('\t')[3]+'-'+line.split('\t')[4]+line.split('\t')[6]
    return dic_out

dic_annotNt=annotfile(Nr_Nt)
dic_annotNb=annotfile(Nr_Nb)
dic_Ntlen=read_fasta.fasta_len(Nta_db)
dic_Nblen=read_fasta.fasta_len(Nbe_db)
dic_Ntloci=gff2loci(Nt_gff)
dic_Nbloci=gff2loci(Nb_gff)

dic_TAIR={}
with open(Atpep,'r') as file_in:
    for line in file_in:
        if re.match(r'>',line):
            geneId=re.match(r'>(\S+)',line).group(1).upper()
            tannot=line.split('|')[2].strip()
            dic_TAIR[geneId]=tannot
dic_slyannot={}
with open(Slpep,'r') as file_in:
    for line in file_in:
        if re.match(r'>',line):
            geneId=re.match(r'>(\S+)',line).group(1)
            transId=re.match(r'Solyc\d+g\d+\.\d+\.\d+',geneId).group()
            tannot=re.sub(r'^\S+\s+','',line.strip())
            tannot=re.sub(r'\(.+','',tannot)
            dic_slyannot[transId]=tannot 

dic_final={}
dic_id2score={}
with open(Blastf,'r') as file_in:
    for line in file_in:
        type1=line.split()[0].split('|')[0]
        type2=line.split()[1].split('|')[0]
        if re.search(r'N',type1+type2) and (re.search(r'ATH\|AT\dG\d+',line) or re.search(r'SLY\|Solyc\d+g\d+',line)):
            tobaccoid=re.search(r'N[tb][ae]\d+g\d+\.\d+',line).group()
            tobaccolen=0
            if re.match(r'Nta',tobaccoid):
                tobaccolen=dic_Ntlen[tobaccoid]
            else:
                tobaccolen=dic_Nblen[tobaccoid]
            matchlen=int(line.split()[3])
            matchratio=float(line.split()[2])
            if (matchlen/tobaccolen>=0.5 and matchratio>=50) or (float(line.strip().split()[-1])>500):
                lineid=re.match(r'\S+\t\S+',line).group()
                anothergeneid=''
                if re.search(r'ATH\|AT\dG',lineid):
                    anothergeneid=re.search(r'AT\dG\d+\.\d+',lineid,flags=re.I).group().upper()
                elif re.search(r'SLY',lineid):
                    anothergeneid=re.search(r'Solyc\d+g\d+\.\d+\.\d+',lineid).group()
                leveltype1=tobaccoid[0:2]+anothergeneid[0:1]
                med=dic_final.setdefault(leveltype1,{})
                med.setdefault(tobaccoid,{})
                leveltype2=tobaccoid+'\t'+anothergeneid
                dic_final[leveltype1][tobaccoid][leveltype2]='Y'
                dic_id2score.setdefault(leveltype2,[])
                dic_id2score[leveltype2].append(line.strip().split()[-1])

dic_id2score2={}
for eachone,eachlist in dic_id2score.items():
    eachlistsort=sorted(eachlist,key=lambda x:-float(x))
    dic_id2score2[eachone]=eachlistsort[0]
out1w=open(out1,'w')
with open(Nta_db,'r') as file_in:
    for line in file_in:
        if re.match(r'>',line):
            transid=re.match(r'>(\S+)',line).group(1)
            athkey=transid[0:2]+'A'
            slykey=transid[0:2]+'S'
            thisannot=''
            athid='-'
            slyid='-'
            if transid in dic_final[athkey].keys():
                athidlist=list(dic_final[athkey][transid].keys())
                athidlistsort=sorted(athidlist,key=lambda x:-float(dic_id2score2[x]))
                athid=athidlistsort[0].split('\t')[1]
                thisannot=dic_TAIR[athid]
            if transid in dic_final[slykey].keys():
                slyidlist=list(dic_final[slykey][transid].keys())
                slyidlistsort=sorted(slyidlist,key=lambda x:-float(dic_id2score2[x]))
                slyid=slyidlistsort[0].split('\t')[1]
                if thisannot=='':
                    thisannot=dic_slyannot[slyid]
            if thisannot=='':
                thisannot=dic_annotNt[transid]
            if thisannot=='' or thisannot=='-':
                thisannot='-'
            out1w.write(transid+'\t'+dic_Ntloci[transid]+'\t'+slyid+'\t'+athid+'\t'+thisannot+'\n')
out1w.close()
out2w=open(out2,'w')
with open(Nbe_db,'r') as file_in:
    for line in file_in:
        if re.match(r'>',line):
            transid=re.match(r'>(\S+)',line).group(1)
            athkey=transid[0:2]+'A'
            slykey=transid[0:2]+'S'
            thisannot=''
            athid='-'
            slyid='-'
            if transid in dic_final[athkey].keys():
                athidlist=list(dic_final[athkey][transid].keys())
                athidlistsort=sorted(athidlist,key=lambda x:-float(dic_id2score2[x]))
                athid=athidlistsort[0].split('\t')[1]
                thisannot=dic_TAIR[athid]
            if transid in dic_final[slykey].keys():
                slyidlist=list(dic_final[slykey][transid].keys())
                slyidlistsort=sorted(slyidlist,key=lambda x:-float(dic_id2score2[x]))
                slyid=slyidlistsort[0].split('\t')[1]
                if thisannot=='':
                    thisannot=dic_slyannot[slyid]
            if thisannot=='':
                thisannot=dic_annotNb[transid]
            if thisannot=='' or thisannot=='-':
                thisannot='-'
            out2w.write(transid+'\t'+dic_Nbloci[transid]+'\t'+slyid+'\t'+athid+'\t'+thisannot+'\n')                
out2w.close()        
          