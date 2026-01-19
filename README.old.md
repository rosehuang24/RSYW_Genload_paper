# Genload

## To prepare variant impact classes using SNPEff and GERP scores
working dir: ```/storage/zhenyingLab/huangruoshi/genload/```

### 1. SNPeff annotation
```
java -Xmx4g -jar $snpeff -c $soft/snpEff/snpEff.config -v GRCg6a.110 \
        /storage/zhenyingLab/huangruoshi/anc_inferrence/116combineVCFs/swapped_aO.${chrm}_${startpos}_${endpos}.vcf.gz \
        > tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out
```
### 2. Extract modifier, synonymous (low) and modifier

```
#although we actually didn't use HIGH category
grep HIGH tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | awk '{print$1"\t"$2-1"\t"$2}' \
      > tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.HIGH.pos.bed

grep MODERATE tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | grep -v HIGH  |awk '{print$1"\t"$2-1"\t"$2}' \
      > tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.MODERATE.pos.bed

grep LOW tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | grep synonymous_variant | grep -v HIGH | grep -v MODERATE  | awk '{print$1"\t"$2-1"\t"$2}' \
      > tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.LOW.synonymous_only.bed

grep MODIFIER tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | grep -v HIGH | grep -v MODERATE | grep -v LOW | awk '{print$1"\t"$2-1"\t"$2}' \
      >  tmp_dir/116.swapped.${chrm}_${startpos}_${endpos}.MODIFIER.pos.bed

```
### 3. Further standardize the impact into deleterious, synonymous and neutral 

#### Deleterious mutations: missense mutations with miyata score >= 1.85

```
python3 $scriptDIR/miyata.del.cutoff.py -I midfiles.miyata.input -O 116.aO_miyata_deleterious.pos.bed -M $parentDIR/Miyata_score/miyata_score_tab.txt 
```

If HIGH are included, then it is included. No miyata process shall be apply to HIGHs because not all of them are protein coding


#### Synonymous mutations:

```
cat tmp_dir/116.swapped.*.LOW.synonymous_only.bed| sort -k 1,1n -k 2,2n > 116.synonymous.pos.bed
```

#### Neutral mutations: 
preperation dir: ```/storage/zhenyingLab/huangruoshi/genload/neutral_region```

Based on 20240116_RH.pptx and 20240122_RH.pptx, we used UCSC_rep coordinates and cds feature from gff3 files. 

It is slightly different from previous pilot run because the release versions of gff3 files are different.

First, we define neutral region:
1. At least 50 kb from CDS region (gff annotation)
2. No three consecutive sites with PC>0.5 with 100 b flanking region
3. no repetitive region (UCSC annotation, same in variant-calling/filtering procedure)
   
```
grep -v "#" $ref/Gallus_gallus_gca000002315v5.GRCg6a.110.gff3  \
|  awk '($3=="CDS") && ($1<29) {print$1"\t"$4"\t"$5}'  \
| sort -k 1,1n -k 2,2n \
| awk '{print$1"\t"$2-50000"\t"$3+50000}' \
> gff3.release110.CDS.50kbflanking.bed

parallel --colsep '\t' -q awk '($1=="{1}") { if($3 > {2}) print $1"\t"$2"\t{2}"; else print $0}'  gff3.release110.CDS.50kbflanking.bed  \
:::: /storage/zhenyingLab/huangruoshi/20211122_nh/genome.from.dna_rm \
| awk '{ if ($2 < 0) print$1"\t0\t"$3; else print$0}' | sort -k 1,1n -k 2,2n | bedtools merge -i - \
> gff3.release110.CDS.50kbflanking.inrange.bed

cat \
$ref/repetitive_auto.bed \
/storage/zhenyingLab/huangruoshi/20211122_nh/PC_larger_than_0.5_length3_100b_flanking.merged.tailed.bed \
gff3.release110.CDS.50kbflanking.inrange.bed \
| sort -k 1,1n -k 2,2n | bedtools merge -i - | awk '($1<29) {print$0}' \
>  dispose.PC_rep_cds_50kbflanking.bed

bedtools complement -i dispose.PC_rep_cds_50kbflanking.bed -g /storage/zhenyingLab/huangruoshi/20211122_nh/genome.from.dna_rm \
| awk '($1<29) {print$0}' \
> neutral.bed

```

All sites in ```neutral.bed``` are considered neutral

