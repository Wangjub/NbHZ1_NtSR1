import re,os,sys,numpy

directory=sys.argv[1]
filelist=sys.argv[2]
#numpy.std(list,ddof=1)
#numpy.mean(list)求均值
for eachunit in filelist.split(';'):
    fileout=eachunit.split(':')[0]
    sampleList=eachunit.split(':')[1].split('#')
    label=sampleList.pop(-1).split(',')
    dic_express={}
    dic_label2sample={}
    for i in range(len(label)):
        dic_label2sample[label[i]]=sampleList[i]
    dic_genelist={}    
    for fileunit in sampleList:
        for eachfile in fileunit.split(','):
            with open(directory+'/'+eachfile+'_express.matrix','r') as file_in:
                file_in.readline()
                for line in file_in:
                    geneId=line.split()[0]
                    dic_genelist[geneId]=1
                    geneex=line.strip().split()[-1]
                    dic_express.setdefault(eachfile,{})
                    dic_express[eachfile][geneId]=geneex
    filewrite=open(fileout+'_express.matrix','w')
    filewrite.write('gene\t'+'\t'.join(label)+'\n')
    for eachgene in sorted(list(dic_genelist.keys())):
        eachline=''
        for eachtype in label:
            label2sample=dic_label2sample[eachtype]
            geneexlist=[]
            for eachone in label2sample.split(','):
                geneexlist.append(float(dic_express[eachone][eachgene]))
            genemean=round(numpy.mean(geneexlist),2)
            genestd=round(numpy.std(geneexlist),2)
            eachline+=str(genemean)+'#'+str(genestd)+'\t'
        filewrite.write(eachgene+'\t'+eachline.strip()+'\n')
    
    filewrite.close()



