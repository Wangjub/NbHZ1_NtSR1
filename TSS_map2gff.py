import re,os,sys

TSSmap=sys.argv[1]
GFFfile=sys.argv[2]
outfile=sys.argv[3]



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
with open(TSSmap,'r') as file_in:
    for line in file_in:
        chro2=line.split()[2]
        strand2=line.split()[1]
        thisloci=''
        if strand2=='+':
            thisloci=line.split()[3]
        elif strand2=='-':
            thisloci=str(int(line.split()[3])+len(line.split()[4]))
        if chro2+'\t'+strand2 in dic_loci.keys() and thisloci in dic_loci[chro2+'\t'+strand2].keys():
            tseq=line.split()[4]
            matchs=int(line.split()[3])+1
            matche=matchs+len(tseq)
            chromo=line.split()[2]
            seqstrand=line.split()[1]
            thiskey=chromo+'\t'+str(matchs)+'\t'+str(matche)+'\t'+seqstrand
            dic_TSS.setdefault(thiskey,{})
            dic_TSS[thiskey][tseq]='Y'
 
number=0            
out=open(outfile,'w')
for eachloci in sorted(list(dic_TSS.keys()),key=lambda x:(x.split()[0] , int(x.split()[1]))):
    chro3=eachloci.split()[0]
    start3=eachloci.split()[1]
    end3=eachloci.split()[2]
    strand3=eachloci.split()[3]
    
    for eachseq in dic_TSS[eachloci].keys():
        number+=1
        out.write(chro3+'\tTSSmap\tTSS\t'+start3+'\t'+end3+'\t.\t'+strand3+'\t.\tID=TSS'+str(number)+'Name=TSS'+str(number)+'\n')
out.close() 