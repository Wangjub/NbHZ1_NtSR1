import re,os,sys
from python_bymyself import read_fasta

Blastf=sys.argv[1]
Nbe_db=sys.argv[2]
NBEfa=sys.argv[3]
output=sys.argv[4]

dic_Nblen=read_fasta.fasta_len(Nbe_db)


dic_final={}
dic_id2score={}
with open(Blastf,'r') as file_in:
    for line in file_in:
        type1=line.split()[0].split('|')[0]
        type2=line.split()[1].split('|')[0]
        if re.search(r'NBHZ',type1+type2) and (re.search(r'NbQld',line) or re.search(r'NBLab',line)):
            tobaccoid=re.search(r'N[tb][ae]\d+g\d+\.\d+',line).group()
            tobaccolen=dic_Nblen[tobaccoid]
            matchlen=int(line.split()[3])
            matchratio=float(line.split()[2])
            if (matchlen/tobaccolen>=0.7 and matchratio>=70) or (float(line.strip().split()[-1])>500):
                lineid=re.match(r'\S+\t\S+',line).group()
                anothergeneid=''
                leveltype1=''
                if re.search(r'NbQld',lineid):
                    anothergeneid=re.search(r'NbQ\d+g\d+\.\d+',lineid,flags=re.I).group()
                    leveltype1=tobaccoid[0:2]+'NbQld'
                elif re.search(r'NBLab',lineid):
                    anothergeneid=re.search(r'NbL\d+g\d+\.\d+',lineid,flags=re.I).group()
                    leveltype1=tobaccoid[0:2]+'NBLab'
                
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
outfile=open(output,'w')
#print(dic_final['NbNB101']['Niben101Scf03546g00006.1'])
with open(NBEfa,'r') as file_in:
    for line in file_in:
        if re.match(r'>',line):
            transid=line.strip().split('|')[-1]
            NbQLDid='-'
            NbLabid='-'
            if transid in dic_final['NbNbQld'].keys():
                qldidlist=list(dic_final['NbNbQld'][transid].keys())
                qldidlistsort=sorted(qldidlist,key=lambda x:-float(dic_id2score2[x]))
                NbQLDid=qldidlistsort[0].split('\t')[1]
            if transid in dic_final['NbNBLab'].keys():
                labidlist=list(dic_final['NbNBLab'][transid].keys())
                labidlistsort=sorted(labidlist,key=lambda x:-float(dic_id2score2[x]))
                NbLabid=labidlistsort[0].split('\t')[1]
            outfile.write(transid+'\t'+NbQLDid+'\t'+NbLabid+'\n')

