# WORKFLOW

## Categories
**Coding**: LoF, missense deleterious, synonymous

SnpEff + miyata score

**Non-coding**: non-coding deleterious, non-coding neutral

GERP score

## 1.1. Coding

use the polarized info
```
java -Xmx4g -jar $snpeff -c $soft/snpEff/snpEff.config -v GRCg6a.110 /storage/zhenyingLab/huangruoshi/anc_inferrence/RSYW.aO.swapped.recode.vcf  > tmp_dir/RSYW.aO.SnpEff_ann.out
```
check for pseudogenes
```
#code for checking (they have same counts)awk :
#wc -l RSYW.aO.MODERATE.noheader.vcf
#awk -F ',' '{print$1}' RSYW.aO.MODERATE.noheader.vcf | grep -v pseudo | wc -l
```

### LOF:
Due to our polarization processing nature, the stop lost and stop gain is opposite.
The actual stop lost is not considered here, therefore the "stop gain" ought to be removed

```
grep LOF tmp_dir/RSYW.aO.SnpEff_ann.out > tmp_dir/LOF.noheader.vcf

cd tmp_dir
awk 'length($9)>3' LOF.noheader.vcf  > original.state.tmp 
awk 'length($9)<3' LOF.noheader.vcf  | grep -v stop_gain > swapped_nostopgain.tmp #because if there is one, there will be on the first anno
cat swapped_nostopgain.tmp original.state.tmp  | sort -k 1,1n -k 2,2n |awk '{print$1"\t"$2-1"\t"$2}' > LOF.final.bed

bedtools intersect -a LOF.final.bed -b $anc/RSYW.aO.swapped.pos.bed > ../LoF.pos.bed
```
```length($9)``` can be an indicater for whether the site was swapped.

### Missense deleterious

```
grep missense RSYW.aO.SnpEff_ann.out | grep -v HIGH > RSYW.aO.MODERATE.noheader.vcf
python3 $scriptDIR/miyata.del.cutoff.py -I RSYW.aO.MODERATE.noheader.vcf  -O missense.miyata.del.bed -M $parentDIR/Miyata_score/miyata_score_tab.txt
cp missense.miyata.del.bed ../deleterious.pos.bed
```

### Synonymous

```
grep synonymous_variant RSYW.aO.SnpEff_ann.out  | grep -v HIGH | grep -v MODERATE  | awk '{print$1"\t"$2-1"\t"$2}' > ../synonymous.pos.bed
```

## 1.2. Non-coding

first, non-coding region

#CDS from gff3 (release110)
ref=/storage/zhenyingLab/huangruoshi/chicken_ref
```
#grep -v "#" $ref/Gallus_gallus_gca000002315v5.GRCg6a.110.gff3  |  awk '($3=="CDS") && ($1<29) {print$1"\t"$4-1"\t"$5}'  | sort -k 1,1n -k 2,2n | bedtools merge -i - > gff3.release110.CDS.bed
bedtools subtract -a gff3.release110.CDS.bed -b $ref/repetitive_auto.bed | sort -k 1,1n -k 2,2n > gff3.release110.CDS.no_rep.bed

grep -v "#" $ref/Gallus_gallus_gca000002315v5.GRCg6a.110.gff3  \
          |  awk '($3=="CDS") && ($1<29) {print$1"\t"$4-1"\t"$5}'  \
          | sort -k 1,1n -k 2,2n | awk '{print$1"\t"$2-50000"\t"$3+50000}' > gff3.release110.CDS.50kbflanking.bed

#make sure it doesn't go beyond chrm boundary
# {1}: chrm. {2}: the size (larger end) of that chrm
parallel --colsep '\t' -q awk '($1=="{1}") { if($3 > {2}) print $1"\t"$2"\t{2}"; else print $0}'  gff3.release110.CDS.50kbflanking.bed  \
          :::: /storage/zhenyingLab/huangruoshi/20211122_nh/genome.from.dna_rm \
          | awk '{ if ($2 < 0) print$1"\t0\t"$3; else print$0}' \
          | sort -k 1,1n -k 2,2n \
          | bedtools merge -i - > gff3.release110.CDS.50kbflanking.inrange.bed

cat $ref/repetitive_auto.bed gff3.release110.CDS.50kbflanking.inrange.bed | sort -k 1,1n -k 2,2n | bedtools merge -i - | awk '($1<29) {print$0}' >  dispose.rep_cds_50kbflanking.bed

bedtools complement -i dispose.rep_cds_50kbflanking.bed  -g /storage/zhenyingLab/huangruoshi/20211122_nh/genome.from.dna_rm | awk '($1<29) {print$0}' > non_CDS_flanking_no_rep.bed
```

The non-coding region: ```non_CDS_flanking_no_rep.bed```


