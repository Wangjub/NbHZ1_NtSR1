import re,os,sys

direct1=sys.argv[1]

for eachfile in os.listdir(direct1):
    if re.search(r'\.id$',eachfile):
        allfile=direct1+'/'+eachfile
        expfile=direct1+'/'+re.sub(r'\.id','',eachfile)+'.matrix'
        allnum=0
        #expnum=0
        up= 0
        down=0
        noch= 0
        idList=[]
        dic_idtype={}
        with open(allfile,'r') as file_in:
            for line in file_in:
                if re.search(r'\S+',line):
                    allnum+=1
                    idList.append(line.strip())
        with open(expfile,'r') as file_in:
            file_in.readline()
            for line in file_in:
                if re.search(r'\S+',line):
                    thisid=re.sub(r'\"','',line.split()[0])
                    express1=float(line.split()[-2])
                    express2=float(line.split()[-1])
                    if express1>=express2:
                        if express2==0:
                            down+=1
                            dic_idtype[thisid]='D'
                        elif express1/express2>=2:
                            down+=1
                            dic_idtype[thisid]='D'
                        else:
                            noch+=1
                            dic_idtype[thisid]='ns'
                    elif express1< express2:
                        if express1==0:
                            up+=1
                            dic_idtype[thisid]='U'
                        elif express2/express1>=2:
                            up+=1
                            dic_idtype[thisid]='U'
                        else:
                            noch+=1
                            dic_idtype[thisid]='ns'
                    
        nonnum=allnum-up-down-noch
        com=open('lin.R','w')
        com.write('pdf(file="'+re.sub(r'\.id','',eachfile)+'.pie.pdf",width=5,height=5)\nlibrary(ggplot2)\n')
        com.write('freq=c('+str(up)+', '+str(down)+', '+str(noch)+', '+str(nonnum)+')\n')
        com.write('mydf <- data.frame(causes=c("up","down","un-ch", "non-ex"),freq=c('+str(up)+', '+str(down)+', '+str(noch)+', '+str(nonnum)+'),share= 100*c('+str(up)+'/sum(freq), '+str(down)+'/sum(freq), '+str(noch)+'/sum(freq), '+str(nonnum)+'/sum(freq)))\n')
        com.write('ggplot(mydf, aes("", share, fill = causes)) +geom_bar(width =1, size =1, color = "white", stat = "identity") + coord_polar("y") +  geom_text(aes(label = paste(freq)), position = position_stack(vjust = 0.5), size=10, col="black") + labs(x = NULL, y = NULL, fill = NULL, title = "") +  guides(fill = guide_legend(reverse = TRUE)) +scale_fill_manual(values = c("green", "#BDBDBD","#6699CC","red")) +theme_classic() + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),legend.text=element_text(size=20))\n')
        com.write('dev.off()\n')
        com.close()
        os.system('Rscript lin.R')
        os.remove('lin.R')
        
        id_ex=open(re.sub(r'\.id','',eachfile)+'.gene.note','w')
        for eachgene in sorted(idList):
            if eachgene in dic_idtype.keys():
                id_ex.write(eachgene+'\t'+dic_idtype[eachgene]+'\n')
            else:
                id_ex.write(eachgene+'\tne\n')
        id_ex.close()

