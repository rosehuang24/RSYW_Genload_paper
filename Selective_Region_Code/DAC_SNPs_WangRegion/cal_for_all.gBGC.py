import argparse
import gzip
parser = argparse.ArgumentParser()
parser.add_argument("-I","--input", help="input file, zipped or unzipped vcf",required=True)
parser.add_argument("-O","--output", help="output file name",required=True)
parser.add_argument("-n","--alleles", help="number of alleles in the samping population.",type=int, required=True)
parser.add_argument("-p","--population", help="name for population, if not provided, will use input name as default",required=False)

args = parser.parse_args()
if args.population:
    name=str(args.population)
else:
    name=str(args.input)

#set number of sites
unbiased_sites=0
SW_sites=0
WS_sites=0
#set number of derived alleles
unbiased_vars=0
SW_vars=0
WS_vars=0

outh = open(args.output, 'w')

unbiased_ls=["AT","TA","GC","CG"]
SW_ls=["GT","GA","CT","CA"]
WS_ls=["AG","TG","AC","TC"]

with (gzip.open if args.input.endswith(".gz") else open)(args.input, "rt") as inh:
    for lines in inh:
        if not lines.startswith("CHROM"):
#            pos=lines.strip().split()[0]+"\t"+str(int(lines.strip().split()[1])-1)+"\t"+lines.strip().split()[1]
            ref=lines.strip().split()[3].split(":")[0]
            alt=lines.strip().split()[4].split(":")[0]
            code=ref+alt
            altfreq=float(lines.strip().split()[4].split(":")[1])
            if code in unbiased_ls:
                unbiased_sites+=1
                unbiased_vars+=args.alleles*altfreq
            elif code in SW_ls:
                SW_sites+=1
                SW_vars+=args.alleles*altfreq
            elif code in WS_ls:
                WS_sites+=1
                WS_vars+=args.alleles*altfreq

outh.write(name+"\tWWSS\t"+str(unbiased_sites)+"\t"+str(unbiased_vars)+"\t"+str(unbiased_vars/unbiased_sites)+"\n")
outh.write(name+"\tSW\t"+str(SW_sites)+"\t"+str(SW_vars)+"\t"+str(SW_vars/SW_sites)+"\n")
outh.write(name+"\tWS\t"+str(WS_sites)+"\t"+str(WS_vars)+"\t"+str(WS_vars/WS_sites)+"\n")

inh.close()
outh.close()
