import re,os,sys
direct1=sys.argv[1]
species=sys.argv[2]


dic_specie2id={}
dic_num={}
dic_orth={}
for eachfile in os.listdir(direct1):
    if re.search(r'^RAxML\_bestTree',eachfile):
        with open(direct1+'/'+eachfile,'r') as file_in:
            thistree=file_in.readline()
            matchlist=re.findall(r'\((Ngla|Nobt|Ntom|Noto|Npan|Nkni|Nund|Nsyl|Natt|Nbet|Ntab|Nrus)\:\d+\.\d+\,(Ngla|Nobt|Ntom|Noto|Npan|Nkni|Nund|Nsyl|Natt|Nbet|Ntab|Nrus)\:\d+\.\d+\)',thistree)
            for eachone in matchlist:
                if species in eachone:
                    species1=eachone[0]
                    species2=eachone[1]
                    thiskey="\t".join(sorted([species1,species2]))
                    dic_orth.setdefault(thiskey,[])
                    orthid=re.sub(r'_tree','',eachfile.split('.')[1])
                    dic_orth[thiskey].append(orthid)
                    if thiskey in dic_num.keys():
                        dic_num[thiskey]+=1
                    else:
                        dic_num[thiskey]=1
    if re.search(r'^OG\d+\.fasta$',eachfile):
        orthfamily=eachfile.split('.')[0]
        with open(direct1+'/'+eachfile,'r') as file_in2:
            for line2 in file_in2:
                if re.match(r'>',line2):
                    med=dic_specie2id.setdefault(orthfamily,{})
                    spe2=re.sub(r'>','',line2.strip()).split('|')[0]
                    spe2id=re.sub(r'>','',line2.strip()).split('|')[1]
                    med.setdefault(spe2,[])
                    dic_specie2id[orthfamily][spe2].append(spe2id)
for eachkey in dic_num.keys():
    print(eachkey,end='\t')
    print(dic_num[eachkey])
    print("\n".join(dic_orth[eachkey]))
    #thisspe=eachkey.split()[1]
    #for eachorth in dic_orth[eachkey]:
    #    print('\t'.join(dic_specie2id[eachorth][species]))
    #print("\n".join(dic_orth[eachkey]))

