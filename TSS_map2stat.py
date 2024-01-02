import re,os,sys

TSSmap=sys.argv[1]
GFFfile=sys.argv[2]
outfile=sys.argv[3]

def statchange(filein,fileout):
    df_test= pd.read_table(filein,sep='\t',index_col=0,header=None)
    df_test.loc["sum"]=df_test.apply(lambda x:sum(x),axis=0)
    df_test['ave']=df_test[1]/df_test.loc['sum'][1]
    df_test.columns=['value','average']
    #df_test.drop([len(df_test)-1],inplace=True)
    df_test = df_test[:-1]
    df_test.to_csv(fileout,sep='\t',index=1)

dic_loci={}
dic_gff={}
with open(GFFfile,'r') as file_in:
    for line in file_in:
        if re.search(r'\tCDS\t',line):
            transId=re.search(r'Parent=(\S+)',line).group(1)
            dic_gff.setdefault(transId,[])
            dic_gff[transId].append(line.split()[0]+'\t'+line.split()[6]+'\t'+line.split()[3]+'\t'+line.split()[4])
for eachid,list2 in dic_gff.items():
    tchro=list2[0].split()[0]
    tstrand=list2[0].split()[1]
    listsort=[]
    CDSstart=0
    CDSend=0
    listsort=sorted(list2,key=lambda x:int(x.split('\t')[2]))
    if tstrand=='+':
        CDSstart=int(listsort[0].split()[2])
        CDSend=listsort[-1].split()[3]
        for i in range(1,2001,1):
            med1=dic_loci.setdefault(tchro+'\t'+tstrand,{})
            med1.setdefault(str(CDSstart-i),[])
            dic_loci[tchro+'\t'+tstrand][str(CDSstart-i)].append(str(i))
    elif tstrand=='-':
        CDSend=listsort[0].split()[2]
        CDSstart=int(listsort[-1].split()[3])
        for i in range(1,2001,1):
            med1=dic_loci.setdefault(tchro+'\t'+tstrand,{})
            med1.setdefault(str(CDSstart+i),[])
            dic_loci[tchro+'\t'+tstrand][str(CDSstart+i)].append(str(i))
     


dic_TSS={}
num=0
with open(TSSmap,'r') as file_in:
    for line in file_in:
        num+=1
        #print(str(num),end='\r')
        chro2=line.split()[2]
        strand2=line.split()[1]
        thisloci=''
        if strand2=='+':
            thisloci=line.split()[3]
        elif strand2=='-':
            thisloci=str(int(line.split()[3])+len(line.split()[4]))
        if chro2+'\t'+strand2 in dic_loci.keys() and thisloci in dic_loci[chro2+'\t'+strand2].keys():
            tread=int(line.split()[0].split('_')[2])
            toTSS=sorted(dic_loci[chro2+'\t'+strand2][thisloci],key=lambda x:-int(x))[0]
            if toTSS in dic_TSS.keys():
                dic_TSS[toTSS]+=tread
            else:
                dic_TSS[toTSS]=tread
            
out=open(outfile,'w')            
for eachnum in sorted(list(dic_TSS.keys()),key=lambda x:int(x)):
    out.write(eachnum+'\t'+str(dic_TSS[eachnum])+'\n')
out.close()    
statchange(outfile,outfile+'.change')   