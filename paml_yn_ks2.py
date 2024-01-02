import re,os,sys

orthcds=sys.argv[1]  #OG11017.fasta
dnds_dir=sys.argv[2]
spec    =sys.argv[3]

dic_unit={}

for eachfile in os.listdir(orthcds):
    if re.search(r'^OG\d+\.fasta$',eachfile):
        orthid=re.sub(r'\.fasta','',eachfile)
        listid=[]
        dic_num2spe={}
        with open(orthcds+'/'+eachfile,'r') as file_in:
            for line in file_in:
                if re.match(r'>',line):
                    specid=re.match(r'>(\S+)',line).group(1).split('|')[0]
                    listid.append(specid)
        idlen=len(listid)
        for i in range(0,idlen-1):
            for j in range(i+1,idlen,1):
                dic_num2spe[str(i+1)+'\t'+str(j+1)]=listid[i]+'\t'+listid[j]
                dic_num2spe[str(j+1)+'\t'+str(i+1)]=listid[i]+'\t'+listid[j]
        with open(dnds_dir+'/'+orthid+'.txt','r') as file_in:
            for line in file_in:
                if re.match(r'\d',line.strip()) and len(line.strip().split())==13:
                    num1=line.strip().split()[0]
                    num2=line.strip().split()[1]
                    tks=line.strip().split()[-3]
                    if re.search(r'nan',tks):
                        pass
                    else:
                        #print(orthid)
                        spekey=dic_num2spe[num1+'\t'+num2]
                        if re.search(spec,spekey):
                            dic_unit.setdefault(spekey,[])
                            dic_unit[spekey].append(tks)
print('species\tKs')
for eachname,eachlist in dic_unit.items():
    for eachnum in eachlist:
        print(re.sub(r'\s+','-',eachname)+'\t'+eachnum)
                

