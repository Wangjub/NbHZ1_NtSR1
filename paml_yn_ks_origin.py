import re,os,sys

orthcds=sys.argv[1]  #OG11017.fasta
dnds_dir=sys.argv[2]
spec    =sys.argv[3]
orth_list=sys.argv[4]

dic_orth={}
with open(orth_list,'r') as file_in:
    spec2=''
    spec2num=''
    for line in file_in:
        if re.match(r'N',line):
            spec2=line.split()[0]+'-'+line.split()[1]
            spec2num=line.strip().split()[-1]
        else:
            dic_orth.setdefault(spec2+'\t'+spec2num,[])
            dic_orth[spec2+'\t'+spec2num].append(line.strip())


dic_orth_final={}
keysort=sorted(list(dic_orth.keys()),key=lambda x:-int(x.split()[-1]))

for eachp in dic_orth[keysort[0]]:
    dic_orth_final[eachp]=keysort[0].split()[0]
for eachp in dic_orth[keysort[1]]:
    dic_orth_final[eachp]=keysort[1].split()[0]

       

dic_unit={}
for eachfile in os.listdir(orthcds):
    if re.search(r'^OG\d+\.fasta$',eachfile):
        orthid=re.sub(r'\.fasta','',eachfile)
        if orthid in dic_orth_final.keys():
            twospec=dic_orth_final[orthid]
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
                            listtrans=spekey.split()
                            if listtrans[0]+'-'+listtrans[1] == twospec or listtrans[1]+'-'+listtrans[0] == twospec:
 #                               print(spekey)
                                if re.search(spec,spekey):
                                    dic_unit.setdefault(spekey,[])
                                    dic_unit[spekey].append(tks)
print('species\tKs')
for eachname,eachlist in dic_unit.items():
    for eachnum in eachlist:
        print(re.sub(r'\s+','-',eachname)+'\t'+eachnum)
                

