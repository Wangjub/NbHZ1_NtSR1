import re,os,sys
from python_bymyself import read_fasta
from collections import Counter
import pandas as pd

blastnf=sys.argv[1]
CDSfile=sys.argv[2]
transfi=sys.argv[3]

dic_CDSlen=read_fasta.fasta_len(CDSfile)
dic_tranlen=read_fasta.fasta_len(transfi)

def statchange(filein,fileout):
    df_test= pd.read_table(filein,sep='\t',index_col=0,header=None)
    df_test.loc["sum"]=df_test.apply(lambda x:sum(x),axis=0)
    df_test['ave']=df_test[1]/df_test.loc['sum'][1]
    df_test.columns=['value','average']
    #df_test.drop([len(df_test)-1],inplace=True)
    df_test = df_test[:-1]
    df_test.to_csv(fileout,sep='\t',index=1)


dic_exists={}
dic_final_PB={}
dic_final_RNA={}
with open(blastnf,'r') as file_in:
    for line in file_in:
        transid=line.split()[0]
        if transid in dic_exists.keys():
            pass
        else:
            if re.match(r'PB',transid):
                matchscore=float(line.split()[2])
                matchlen=int(line.split()[3])
                matchgap=int(line.split()[4])
                matchmis=int(line.split()[5])
                ss=line.split()[8]
                se=line.split()[9]
                if matchscore>95 and matchlen>300 and matchgap<=3 and matchmis<=3 and int(ss)<int(se):
                    qs=line.split()[6]
                    qe=line.split()[7]
                    thistransid=line.split()[1]
                    if int(ss)==1 and int(qs)>1:
                        med=dic_final_PB.setdefault('TSS',{})
                        geneId=re.sub(r'\.\d+$','',thistransid)
                        med.setdefault(geneId,[])
                        dic_final_PB['TSS'][geneId].append(int(qs)-1)
                    elif int(se)==dic_CDSlen[thistransid] and int(qe)>int(se):
                        TTSlen=int(qe)-int(se)-1
                        med1=dic_final_PB.setdefault('TTS',{})
                        geneId=re.sub(r'\.\d+$','',thistransid)
                        med1.setdefault(geneId,[])
                        dic_final_PB['TTS'][geneId].append(TTSlen)
            elif re.match(r'MSTRG',transid):
                matchscore=float(line.split()[2])
                matchlen=int(line.split()[3])
                matchgap=int(line.split()[4])
                matchmis=int(line.split()[5])
                ss=line.split()[8]
                se=line.split()[9]
                if matchscore>95 and matchlen>300 and matchgap<=3 and matchmis<=3 and int(ss)<int(se):
                    qs=line.split()[6]
                    qe=line.split()[7]
                    thistransid=line.split()[1]
                    if int(ss)==1 and int(qs)>1:
                        med=dic_final_RNA.setdefault('TSS',{})
                        geneId=re.sub(r'\.\d+$','',thistransid)
                        med.setdefault(geneId,[])
                        dic_final_RNA['TSS'][geneId].append(int(qs)-1)
                    elif int(se)==dic_CDSlen[thistransid] and int(qe)>int(se):
                        TTSlen=int(qe)-int(se)-1
                        med1=dic_final_RNA.setdefault('TTS',{})
                        geneId=re.sub(r'\.\d+$','',thistransid)
                        med1.setdefault(geneId,[])
                        dic_final_RNA['TTS'][geneId].append(TTSlen)
        dic_exists[transid]='Y'

print('PacBio TSS gene number: ',end='')
print(len(dic_final_PB['TSS'].keys()))
print('PacBio TTS gene number: ',end='')
print(len(dic_final_PB['TTS'].keys()))
print('RNA-seq TSS gene number: ',end='')
print(len(dic_final_RNA['TSS'].keys()))
print('RNA-seq TTS gene number: ',end='')
print(len(dic_final_RNA['TTS'].keys()))

fileout1=open('PBTSS.stat','w')
PBlistTSS=[]
for sublist in dic_final_PB['TSS'].values():
    for eachone in sublist:
        PBlistTSS.append(eachone)
PBlistnum=Counter(PBlistTSS)
for i in range(1,max(PBlistTSS)+1):
    if i in PBlistnum.keys():
        fileout1.write(str(i)+'\t'+str(PBlistnum[i])+'\n')
    else:
        fileout1.write(str(i)+'\t0\n')
fileout1.close()
statchange('PBTSS.stat','PBTSS.stat.change')
id1=open('PB_TSS.id','w')
id1.write('\n'.join(list(dic_final_PB['TSS'].keys())))
id1.close()


fileout2=open('PBTTS.stat','w')
PBlistTTS=[]
for sublist in dic_final_PB['TTS'].values():
    for eachone in sublist:
        PBlistTTS.append(eachone)
PBlistnum2=Counter(PBlistTTS)
for i in range(1,max(PBlistTTS)+1):
    if i in PBlistnum2.keys():
        fileout2.write(str(i)+'\t'+str(PBlistnum2[i])+'\n')
    else:
        fileout2.write(str(i)+'\t0\n')
fileout2.close()
statchange('PBTTS.stat','PBTTS.stat.change')
id2=open('PB_TTS.id','w')
id2.write('\n'.join(list(dic_final_PB['TTS'].keys())))
id2.close()

fileout3=open('RNATSS.stat','w')
RNAlistTSS=[]
for sublist in dic_final_RNA['TSS'].values():
    for eachone in sublist:
        RNAlistTSS.append(eachone)
PBlistnum3=Counter(RNAlistTSS)
for i in range(1,max(RNAlistTSS)+1):
    if i in PBlistnum3.keys():
        fileout3.write(str(i)+'\t'+str(PBlistnum3[i])+'\n')
    else:
        fileout3.write(str(i)+'\t0\n')
fileout3.close()
statchange('RNATSS.stat','RNATSS.stat.change')
id3=open('RNA_TSS.id','w')
id3.write('\n'.join(list(dic_final_RNA['TSS'].keys())))
id3.close()

fileout4=open('RNATTS.stat','w')        
RNAlistTTS=[]
for sublist in dic_final_RNA['TTS'].values():
    for eachone in sublist:
        RNAlistTTS.append(eachone)
PBlistnum4=Counter(RNAlistTTS)
for i in range(1,max(RNAlistTTS)+1):
    if i in PBlistnum4.keys():
        fileout4.write(str(i)+'\t'+str(PBlistnum4[i])+'\n')
    else:
        fileout4.write(str(i)+'\t0\n')
fileout4.close()
statchange('RNATTS.stat','RNATTS.stat.change')
id4=open('RNA_TTS.id','w')
id4.write('\n'.join(list(dic_final_RNA['TTS'].keys())))
id4.close()