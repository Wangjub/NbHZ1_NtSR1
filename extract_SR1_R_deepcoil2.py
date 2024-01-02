import re,os,sys

interproscan  =sys.argv[1]
NBS_fa        =sys.argv[2]
TIR_fa        =sys.argv[3]
deepcoil_dir  =sys.argv[4] 


def fasta2list(file):
    outlist=[]
    with open(file,'r') as file_in:
        for line in file_in:
            if re.match(r'>',line):
                geneId=re.sub(r'\.\d+$','',re.match(r'>(\S+)',line).group(1))
                outlist.append(geneId)
    return outlist
    

dic_inter={}

with open(interproscan,'r') as interpro_open:
    for line in interpro_open:
        if re.search(r'Leucine',line,flags=re.I):
            gene_id=re.sub(r'\.\d+$','',line.split()[0])
            dic_inter.setdefault('LRR',{})
            dic_inter['LRR'][gene_id]='Y'

dic_final={}

dic_final['LRR']=list(dic_inter['LRR'].keys())

dic_final['NBS']=fasta2list(NBS_fa)
dic_final['TIR']=fasta2list(TIR_fa)

for eachfile in os.listdir(deepcoil_dir):
    if re.search(r'\.out$',eachfile):
        tlist=[]
        with open(deepcoil_dir+'/'+eachfile,'r') as file_in:
            number=0
            for line in file_in:
                if not re.match(r'aa',line):
                    if float(line.split()[1])>=0.1:
                        number+=1
                    else:
                        tlist.append(number)
                        number=0
        if max(tlist)>=14:
            geneId=re.sub(r'\.out$','',eachfile)
            dic_final.setdefault('CC',[])
            dic_final['CC'].append(geneId)


CTNLlist=set(dic_final['CC']) & set(dic_final['TIR']) & set(dic_final['NBS']) & set(dic_final['LRR'])
TNLlist =set(dic_final['TIR']) & set(dic_final['NBS']) & set(dic_final['LRR']) - set(dic_final['CC'])
CNLlist =set(dic_final['CC']) & set(dic_final['NBS']) & set(dic_final['LRR']) - set(dic_final['TIR'])
CNlist  =set(dic_final['CC']) & set(dic_final['NBS']) - set(dic_final['LRR']) - set(dic_final['TIR'])
NLlist  =set(dic_final['LRR']) & set(dic_final['NBS']) - set(dic_final['CC']) - set(dic_final['TIR'])
TNlist  =set(dic_final['TIR']) & set(dic_final['NBS']) - set(dic_final['CC']) - set(dic_final['LRR'])
Tlist   =set(dic_final['TIR']) - set(dic_final['NBS']) - set(dic_final['CC']) - set(dic_final['LRR'])

if CTNLlist:
    for eachp in sorted(CTNLlist):
        print('CTNL\t'+eachp)
if TNLlist:
    for eachp in sorted(TNLlist):
        print('TNL\t'+eachp)
if CNLlist:
    for eachp in sorted(CNLlist):
        print('CNL\t'+eachp)
if CNlist:  
    for eachp in sorted(CNlist):
        print('CN\t'+eachp)
if NLlist: 
    for eachp in sorted(NLlist):
        print('NL\t'+eachp)
    
if TNlist:
    for eachp in sorted(TNlist):
        print('TN\t'+eachp)

if Tlist:
    for eachp in sorted(Tlist):
        print('T\t'+eachp)

