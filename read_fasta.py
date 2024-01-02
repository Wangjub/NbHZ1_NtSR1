import re
def read_fasta(file):
    protein_dic={}
    with open(file,'r') as protein_open:
        for line in protein_open:
            if line.startswith('>'):
                gene_id=re.sub(r'>','',line).split()[0]
                protein_dic[gene_id]=''
            else:
                protein_dic[gene_id]+=line.strip()
    return protein_dic
    
def read_fasta_gene(file):
    protein_dic={}
    with open(file,'r') as protein_open:
        for line in protein_open:
            if line.startswith('>'):
                gene_id=re.sub(r'>','',line).split()[0]
                gene_id=re.sub(r'\.\d+$','',gene_id)
                protein_dic[gene_id]=''
            else:
                protein_dic[gene_id]+=line.strip()
    return protein_dic
    
def fasta_id_list(file):
    list2=[]
    with open(file,'r') as file_open:
        for line2 in file_open:
            gene_id2=re.sub(r'>','',line2).split()[0]
            if gene_id2 in list2:
                pass
            else:
                list2.append(gene_id2)
    return list2
    
    
def list_2_unique(list_receive):
    list2=[]
    for each in list_receive:
        if each in list2:
            pass
        else:
            list2.append(each)
    return list2

def id_2_annot(file):
    dic_annot={}
    with open(file,'r') as file_open:
        for line2 in file_open:
            if line2.startswith('>'):
                gene_id=re.sub(r'>','',line2).split()[0]
                annot=re.sub(r'>','',line2.strip())
                dic_annot[gene_id]=annot
    return dic_annot  

def fasta_len(file):
    dic_len={}
    with open(file,'r') as file_open:
        for line in file_open:
            if re.match(r'>',line):
                gene_id=re.search(r'>(\S+)',line).group(1)
                dic_len[gene_id]=0
            else:
                dic_len[gene_id]+=len(line.strip())
    return dic_len

def fasta_len_gene(file):
    dic_len={}
    with open(file,'r') as file_open:
        for line in file_open:
            if re.match(r'>',line):
                gene_id=re.search(r'>(\S+)',line).group(1)
                gene_id=re.sub(r'\.\S+$','',gene_id)
                dic_len[gene_id]=0
            else:
                dic_len[gene_id]+=len(line.strip())
    return dic_len