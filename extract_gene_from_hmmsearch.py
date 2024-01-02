import sys,re

protein_db =sys.argv[1]
hmmsearch  =sys.argv[2]
output_file=sys.argv[3]

Fbox_list=[]

hmmsearch_open=open(hmmsearch,'r') 
hmmall_line=hmmsearch_open.read()
hmmsearch_open.close()
hmm_Fbox=re.search(r'(Scores[\s\S]+?\n\n\n)',hmmall_line).group()
for each in hmm_Fbox.split('\n'):
    if re.match(r'\s+\d+',each):
        gene_id=each.strip().split()[8]
        Fbox_list.append(gene_id)

'''
dic_protein={}
with open(protein_db,'r') as protein_open:
    for line in protein_open:
        if re.search(r'>',line): #line.startswith('>')
            gene_id=re.sub(r'>','',line.split('|')[0]) #re.sub(r'old','new',string)
            dic_protein[gene_id]=''
        else:
           dic_protein[gene_id]+=line

output_open=open(output_file,'w')            
for word in Fbox_list:
    output_open.write('>'+word+"\n"+dic_protein[word])
output_open.close()

'''

output_open=open(output_file,'w')
with open(protein_db,'r') as protein_open:
    check=0
    for line in protein_open:
        if re.search(r'>',line):
            gene_id=re.sub(r'>','',line.split()[0])
            if gene_id in Fbox_list:
                output_open.write(line)
                check=1
            else:
                check=0
        elif check==1:
           line2=re.sub(r'\*','',line)
           output_open.write(line2)  
output_open.close()
