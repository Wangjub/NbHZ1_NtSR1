import re,os,sys
direct1=sys.argv[1]
species=sys.argv[2]

orthList=[]
for eachfile in os.listdir(direct1):
    if re.search(r'^RAxML\_bestTree',eachfile):
        with open(direct1+'/'+eachfile,'r') as file_in:
            thistree=file_in.readline()
            if re.search(species,thistree):
                orthList.append(re.search(r'.+?(OG\d+)',eachfile).group(1))

alltype=["Natt","Ngla","Nkni","Nobt","Noto","Npan","Nsyl","Ntom","Nund"]
for eachone in sorted(orthList):
    thisfile=direct1+'/'+eachone+'.fasta'
    dic_final={}
    with open(thisfile,'r') as file_in:
        for line in file_in:
            if re.match(r'>',line):
                dic_final[re.sub(r'>','',line.split('|')[0])]=line.strip().split('|')[-1]

    outp=eachone+'\t'+dic_final[species]+'\t'
    for eachtype in alltype:
        if eachtype in dic_final.keys():
            outp+=dic_final[eachtype]+'\t'
        else:
            outp+='\t'
    outp=re.sub(r'\t$','',outp)
    print(outp)
