import re,os,sys

geneList=sys.argv[1]
nbHZgff=sys.argv[2]
othergff=sys.argv[3]

dic_idlist={}
with open(geneList,'r') as file_in:
    for line in file_in:
        if re.match(r'\S+',line):
            lineNew=re.sub(r'\:\d+','',line.strip())
            list3=lineNew.split()
            geneId=list3.pop(0)
            dic_idlist[geneId]=list3
            
dic_gffnb={}
with open(nbHZgff,'r') as file_in:
    for line in file_in:
        if re.search(r'\tCDS\t',line):
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            med=dic_gffnb.setdefault(transId,{})
            med.setdefault('CDS',[])
            dic_gffnb[transId]['CDS'].append(line.split()[3]+'\t'+line.split()[4])

        if re.search(r'UTR\t',line):
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            med=dic_gffnb.setdefault(transId,{})
            med.setdefault('UTR',[])
            dic_gffnb[transId]['UTR'].append(line.split()[3]+'\t'+line.split()[4])
        if re.search(r'\texon\t',line):
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            med=dic_gffnb.setdefault(transId,{})
            med.setdefault('exon',[])
            dic_gffnb[transId]['exon'].append(line.split()[3]+'\t'+line.split()[4])
        if re.search(r'\tfive_prime_UTR\t',line) and line.split()[6]=='+':
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            med=dic_gffnb.setdefault(transId,{})
            med.setdefault('UTR5',[])
            dic_gffnb[transId]['UTR5'].append(line.split()[3]+'\t'+line.split()[4])
        if re.search(r'\tthree_prime_UTR\t',line) and line.split()[6]=='-':
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            med=dic_gffnb.setdefault(transId,{})
            med.setdefault('UTR5',[])
            dic_gffnb[transId]['UTR5'].append(line.split()[3]+'\t'+line.split()[4])
dic_gffother={}
with open(othergff,'r') as file_in:
    for line in file_in:
        if re.search(r'\tCDS\t',line):
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            transId=re.sub(r'\;.+','',transId)
            med=dic_gffother.setdefault(transId,{})
            med.setdefault('CDS',[])
            dic_gffother[transId]['CDS'].append(line.split()[3]+'\t'+line.split()[4])
        if re.search(r'UTR\t',line):
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            transId=re.sub(r'\;.+','',transId)
            med=dic_gffother.setdefault(transId,{})
            med.setdefault('UTR',[])
            dic_gffother[transId]['UTR'].append(line.split()[3]+'\t'+line.split()[4])
        if re.search(r'\texon\t',line):
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            transId=re.sub(r'\;.+','',transId)
            med=dic_gffother.setdefault(transId,{})
            med.setdefault('exon',[])
            dic_gffother[transId]['exon'].append(line.split()[3]+'\t'+line.split()[4])
        if re.search(r'\tfive_prime_UTR\t',line) and line.split()[6]=='+':
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            transId=re.sub(r'\;.+','',transId)
            med=dic_gffother.setdefault(transId,{})
            med.setdefault('UTR5',[])
            dic_gffother[transId]['UTR5'].append(line.split()[3]+'\t'+line.split()[4])
        if re.search(r'\tthree_prime_UTR\t',line) and line.split()[6]=='-':
            transId=re.search(r'Parent=(\S+)',line.strip()).group(1)
            transId=re.sub(r'\;.+','',transId)
            med=dic_gffother.setdefault(transId,{})
            med.setdefault('UTR5',[])
            dic_gffother[transId]['UTR5'].append(line.split()[3]+'\t'+line.split()[4])

IBSheader='''<?xml version="1.0" encoding="UTF-8"?>
<project version="IBS 1.0" mode="protein" width="800" height="600" rate="1.0">'''
IBSfooter='''</project>'''
def extract_cds(vert,start,end):
    thismodule='''<protein vertical="%s" horizontal="148" locked="true">
    <dash>1</dash>
    <id>CDS</id>
    <start direction="Hide">%s</start>
    <end direction="Hide">%s</end>
    <height size="12" color="255,102,102" colorline="255,102,102">25</height>
    <dash>1</dash>
    <gradient>Horizontal</gradient>
    <texture>00000_0.jpg</texture>
    <textfont />
    </protein>''' % (vert,start,end)
    return thismodule
def extract_utr(vert,start,end):
    thismodule='''<protein vertical="%s" horizontal="148" locked="true">
    <dash>1</dash>
    <id>UTR</id>
    <start direction="Hide">%s</start>
    <end direction="Hide">%s</end>
    <height size="12" color="0,0,255" colorline="0,0,255">25</height>
    <dash>1</dash>
    <gradient>Horizontal</gradient>
    <texture>00000_0.jpg</texture>
    <textfont />
  </protein>''' % (vert,start,end)
    return thismodule
  
def extract_intron(vert,start,end):
    thismodule='''<protein vertical="%s" horizontal="148" locked="true">
    <dash>1</dash>
    <id>intron</id>
    <start direction="Hide">%s</start>
    <end direction="Hide">%s</end>
    <height size="12" color="0,0,0" colorline="255,51,51">15</height>
    <dash>1</dash>
    <gradient>Horizontal</gradient>
    <texture>00000_0.jpg</texture>
    <textfont />
  </protein>''' % (vert,start,end)
    return thismodule
  
def extract_text(deltax,deltay,transId):
    thismodule='''<component type="note">
      <text deltax="%s" deltay="%s" size="25" style="Bold" angle="0" linesize="0" color="0,0,0" linecolor="0,0,0" linedash="false">%s</text>
      <id>%s</id>
      <textfont>Arial</textfont>
    </component>''' % (deltax,deltay,transId,transId)
    return thismodule
def get_intron(exonList):
    exonListSort=sorted(exonList,key=lambda x:int(x.split()[0]))
    intronList=[]
    if len(exonListSort)==1:
        return intronList
    else:
        for i in range(0,len(exonListSort)-1,1):
            introns=int(exonListSort[i].split()[1])+1
            introne=int(exonListSort[i+1].split()[-1])-1
            intronList.append(str(introns)+'\t'+str(introne))
        return intronList
def get_max_min(thisList):
    array1Sort=sorted(thisList,key=lambda x:int(x.split()[0]))
    return array1Sort[0].split()[0],array1Sort[-1].split()[-1]
def retrive_type(nbList,otherList):
    if nbList and otherList:
        nbListSort=sorted(nbList,key=lambda x:int(x.split()[0]))
        otherListSort=sorted(otherList,key=lambda x:int(x.split()[0]))
        nbutrfirst=int(nbListSort[0].split()[1])-int(nbListSort[0].split()[0])
        otherutrfirst=int(otherListSort[0].split()[1])-int(otherListSort[0].split()[0])
        ttype=''
        tlen=0
        if nbutrfirst>=otherutrfirst:
            ttype='other'
            tlen=nbutrfirst-otherutrfirst
        elif nbutrfirst<otherutrfirst:
            ttype='nb'
            tlen=otherutrfirst-nbutrfirst
        return ttype,tlen
    elif nbList:
        nbListSort=sorted(nbList,key=lambda x:int(x.split()[0]))
        nbutrfirst=int(nbListSort[0].split()[1])-int(nbListSort[0].split()[0])
        return 'other',nbutrfirst
    elif otherList:
        otherListSort=sorted(otherList,key=lambda x:int(x.split()[0]))
        otherutrfirst=int(otherListSort[0].split()[1])-int(otherListSort[0].split()[0])
        return 'nb',otherutrfirst
    else:
        return '',''
        
    
for eachid,otherIdList in dic_idlist.items():
    eachidout=open(eachid+'.IBS.xml','w')
    eachidout.write(IBSheader+'\n')
    nbHZintron=get_intron(dic_gffnb[eachid]['exon'])
    nbCDSList=dic_gffnb[eachid]['CDS']
    nbUTRList=[]
    nbUTR5List=[]
    if 'UTR5' in dic_gffnb[eachid].keys():
        nbUTR5List=dic_gffnb[eachid]['UTR5']
    if 'UTR' in dic_gffnb[eachid].keys():
        nbUTRList=dic_gffnb[eachid]['UTR']
    otherUTRlist=[]
    otherfirstId=otherIdList[0]
    if 'UTR5' in dic_gffother[otherfirstId].keys():
        otherUTRlist=dic_gffother[otherfirstId]['UTR5']
    (nb_other,addlen)=retrive_type(nbUTR5List,otherUTRlist)
    (nbmin,nbmax)=get_max_min(nbHZintron+nbCDSList+nbUTRList)
    nbtextdeltax=int((int(nbmax)-int(nbmin))/2)
    nbtextdeltay='44'
    nbIBSList=[]
    for eachintron in nbHZintron:
        eachstart=int(eachintron.split()[0])-int(nbmin)+1
        eachend=int(eachintron.split()[-1])-int(nbmin)+1
        if nb_other=='nb':
            eachstart=eachstart+addlen
            eachend=eachend+addlen
        intronmod=extract_intron('68',str(eachstart),str(eachend))
        nbIBSList.append(intronmod)
    for eachcds in nbCDSList:
        eachstart=int(eachcds.split()[0])-int(nbmin)+1
        eachend=int(eachcds.split()[-1])-int(nbmin)+1
        if nb_other=='nb':
            eachstart=eachstart+addlen
            eachend=eachend+addlen
        cdsmod=extract_cds('68',str(eachstart),str(eachend))
        nbIBSList.append(cdsmod)
    for eachutr in nbUTRList:
        eachstart=int(eachutr.split()[0])-int(nbmin)+1
        eachend=int(eachutr.split()[-1])-int(nbmin)+1
        if nb_other=='nb':
            eachstart=eachstart+addlen
            eachend=eachend+addlen
        utrmod=extract_utr('68',str(eachstart),str(eachend))
        if utrmod:
            nbIBSList.append(utrmod)
    
    lastnbIBS=nbIBSList.pop(-1)
    #print(nbIBSList)
    lastnbList=lastnbIBS.split('\n')
    lastnbone=lastnbList.pop(-1)
    nbtextmo=extract_text(str(nbtextdeltax),nbtextdeltay,eachid)
    nbIBSout='\n'.join(nbIBSList)+'\n'+'\n'.join(lastnbList)+'\n'+nbtextmo+'\n'+lastnbone
    eachidout.write(nbIBSout+'\n')
    
    otherallarray=[]
    for eachotherid in otherIdList:
        otherallarray+=dic_gffother[eachotherid]['CDS']
        if 'UTR' in dic_gffother[eachotherid].keys():
            otherallarray+=dic_gffother[eachotherid]['UTR']
    (othermin,othermax)=get_max_min(otherallarray) 
    for eachotherid in otherIdList:
        otherintronlist=get_intron(dic_gffother[eachotherid]['exon'])
        othercdslist=dic_gffother[eachotherid]['CDS']
        otherutrlist=[]
        if 'UTR' in dic_gffother[eachotherid].keys():
            otherutrlist=dic_gffother[eachotherid]['UTR']
        (eachonemin,eachonemax)=get_max_min(othercdslist+otherutrlist)
        othertextdeltax=int((int(eachonemax)-int(eachonemin))/2)+int(eachonemin)-int(othermin)
        othertextdeltay='28'
        otherIBSList=[]
        for eachintron in otherintronlist:
            eachstart=int(eachintron.split()[0])-int(othermin)+1
            eachend=int(eachintron.split()[-1])-int(othermin)+1
            if nb_other=='other':
                eachstart=eachstart+addlen
                eachend=eachend+addlen
            intronmod=extract_intron('2',str(eachstart),str(eachend))
            otherIBSList.append(intronmod)
        for eachothercds in othercdslist:
            eachstart=int(eachothercds.split()[0])-int(othermin)+1
            eachend=int(eachothercds.split()[-1])-int(othermin)+1
            if nb_other=='other':
                eachstart=eachstart+addlen
                eachend=eachend+addlen
            cdsmod=extract_cds('2',str(eachstart),str(eachend))
            otherIBSList.append(cdsmod)
        for eachutr in otherutrlist:
            eachstart=int(eachutr.split()[0])-int(othermin)+1
            eachend=int(eachutr.split()[-1])-int(othermin)+1
            if nb_other=='other':
                eachstart=eachstart+addlen
                eachend=eachend+addlen
            utrmod=extract_utr('2',str(eachstart),str(eachend))
            otherIBSList.append(utrmod)
        
        lastotherIBS=otherIBSList.pop(-1)
        lastotherList=lastotherIBS.split('\n')
        lastotherone=lastotherList.pop(-1)
        othertextmo=extract_text(str(othertextdeltax),othertextdeltay,eachotherid)
        otherIBSout='\n'.join(otherIBSList)+'\n'+'\n'.join(lastotherList)+'\n'+othertextmo+'\n'+lastotherone
        eachidout.write(otherIBSout+'\n')
    eachidout.write(IBSfooter)