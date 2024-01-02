import re,os,sys
dic_all={}
with open(sys.argv[1],'r') as file_in:
    for line in file_in:
        geneId=re.sub(r'\.t\d+','',line.split()[0])
        types=line.split()[1]
        types=re.sub(r'CTNL','TNL',types)
        dic_all.setdefault(types,{})
        dic_all[types][geneId]=1
idlist=['CLK','CN','CNL','T','TNL','TN','NL','RLP','RLK']
for eachtype in idlist:
    eachout='Niben_'+eachtype+'.id'
    eachwrite=open(eachout,'w')
    eachwrite.write("\n".join(list(dic_all[eachtype].keys())))
    eachwrite.close()
    print(eachtype+'\t',end='')
    print(len(dic_all[eachtype].keys()))
