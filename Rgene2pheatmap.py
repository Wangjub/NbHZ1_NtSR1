import re,os,sys
from python_bymyself import read_fasta

directory1=sys.argv[1]
cuffdiff=sys.argv[2]
R_motif=sys.argv[3]
gene5UTR=sys.argv[4]

dic_genelen=read_fasta.fasta_len_gene(gene5UTR)
dic_Rmotif={}
for eachfile in R_motif.split(','):
    with open(eachfile,'r') as file_in:
        for line in file_in:
            dic_Rmotif[line.strip()]=1

for eachfile in os.listdir(directory1):
    if re.search(r'\.id$',eachfile):
        pdfout=re.sub(r'\.id','',eachfile)+'.pdf'
        os.system('Rscript D:\R_script\geneId2pheatmap2.r '+eachfile+' '+cuffdiff+' cluster_gene:yes;cluster_sample:no;maphigh:red;maplow:blue;cellwidth:20;cellHeight:10;fontsize_row:4;fontsize_col:11;fontsize:20;scale:row;border_color:no;main_name:"Niben CNL pheatmap" '+pdfout)
        #print('Rscript D:\R_script\geneId2pheatmap.r '+eachfile+' '+cuffdiff+' cluster_gene:yes;cluster_sample:no;maphigh:red;maplow:blue;cellwidth:20;cellHeight:10;fontsize_row:4;fontsize_col:11;fontsize:20;scale:row;border_color:no;main_name:"Niben CNL pheatmap" '+pdfout)
        matrixfile=re.sub(r'\.id','',eachfile)+'.normalize.order.matrix'
        Rmotifout=re.sub(r'\.id','',eachfile)+'.normalize.order.matrix.Rmotif'
        Rmotifwrite=open(Rmotifout,'w')
        Rmotifwrite.write("\"R_motif\"\n")
        with open(matrixfile,'r') as file_in:
            file_in.readline()
            for line in file_in:
                geneid=re.sub(r'\"','',line.split()[0])
                Ryesno='0'
                if geneid in dic_Rmotif.keys():
                    Ryesno='1'
                elif geneid in dic_genelen.keys() and dic_genelen[geneid]>15:
                    Ryesno='-1'
                Rmotifwrite.write('"'+geneid+'"\t'+Ryesno+'\n')
                