#input should be SNPeff anno file

import argparse
import gzip
parser = argparse.ArgumentParser()
parser.add_argument("-I","--input", help="input file, zipped or unzipped vcf",required=True)
parser.add_argument("-O","--output", help="output file name",required=True)
parser.add_argument("-M","--miyata", help="file containing mitaya score for all AA changes",required=True)
args = parser.parse_args()

outh = open(args.output, 'w')
mh = open(args.miyata, 'r')
miyata={}
for aapairs in mh:
    aapair="_".join(aapairs.strip().split()[0:2])
    score=float(aapairs.strip().split()[2])
    if score>1.85:
        miyata[aapair]=score

#print(miyata)

with (gzip.open if args.input.endswith(".gz") else open)(args.input, "rt") as inh:
    for lines in inh:
        if not lines.startswith("#"):
            ls=[]
            chrm=lines.strip().split()[0]
            pos=lines.strip().split()[1]
            for field in lines.strip().split()[7].split(","):
                info=str(field).split("|")
                if info[1]=="missense_variant":
                    a1=info[10][2:5]
                    a2=info[10][-3:]
                    aa=(a1+"_"+a2).upper()
                    ls.append(aa) #in case they have different aa change based on different transcript version, although highly unlikely
            #print(chrm+"\t"+pos)
            #print(set(ls))
            uniq_list=set(ls)
            #print(uniq_list)
            for i in uniq_list:
                if i in miyata:
                    outh.write(chrm+"\t"+str(int(pos)-1)+"\t"+pos+"\n")
                    break

inh.close()
outh.close()
mh.close()
