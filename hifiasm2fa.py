import re,os,sys

input_file=sys.argv[1]
output_file=sys.argv[2]

file_out=open(output_file,'w')
with open(input_file,'r') as file_in:
    for line in file_in:
        if re.match(r'S',line):
            list2=line.strip().split()
            file_out.write('>'+list2[1]+' '+list2[3]+' '+list2[4]+'\n'+list2[2]+'\n')
            
file_out.close()

