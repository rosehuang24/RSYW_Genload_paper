# Genload

## To prepare variant impact classes using SNPEff and GERP scores
working dir: ```/storage/zhenyingLab/huangruoshi/genload/```

### 1. SNPeff annotation
```
java -Xmx4g -jar $snpeff -c $soft/snpEff/snpEff.config -v GRCg6a.99 \
        /storage/zhenyingLab/huangruoshi/anc_inferrence/116combineVCFs/swapped.${chrm}_${startpos}_${endpos}.vcf.gz \
        > vcfs/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out
```
### 2. Extract modifier, synonymous (low) and modifier

```
#although we actually didn't use HIGH category
#grep HIGH vcfs/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | awk '{print$1"\t"$2-1"\t"$2}' \
      > snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.HIGH.pos.bed

grep MODERATE vcfs/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | grep -v HIGH  |awk '{print$1"\t"$2-1"\t"$2}' \
      > snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.MODERATE.pos.bed

grep LOW vcfs/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | grep synonymous_variant | grep -v HIGH | grep -v MODERATE  | awk '{print$1"\t"$2-1"\t"$2}' \
      > snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.LOW.synonymous_only.bed

grep MODIFIER vcfs/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | grep -v HIGH | grep -v MODERATE | grep -v LOW | awk '{print$1"\t"$2-1"\t"$2}' \
      >  snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.MODIFIER.pos.bed

```
### 3. Further standardize the impact into deleterious, synonymous and neutral (GERP score and neutral region from ```NH_ZF_CLBL_etc/G-PHOCS_freeman_NH.md```
GERP score: https://ftp.ensembl.org/pub/release-104/compara/conservation_scores/27_sauropsids.gerp_conservation_score/README

#### Deleterious mutations: missense mutations with gerp score > 1
```
parallel bedtools intersect \
        -a snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.MODERATE.pos.bed \
        -b $REFDIR/GERP_score_bins/{}.bed -wa -wb \
        \> GERP_dis/midifles/116.${chrm}_${startpos}_${endpos}.MODERATE.gerp_{} \
        :::: GERPFILES.list
#an example of output file names is 116.1_20000001_40000000.MODERATE.gerp_gerp_scores_-0.6_to_-0.4

cat GERP_dis/midfiles/116.*.MODERATE.gerp_gerp_scores* | awk '($7>=1) {print$1"\t"$2"\t"$3}' | sort -k 1,1n -k 2,2n > 116.deleterious.pos.bed
```
#### Synonymous mutations:

```
cat snpeff_anno/116.swapped.*.LOW.synonymous_only.bed | sort -k 1,1n -k 2,2n > 116.synonymous.pos.bed
```

#### Neutral mutations: 

```
parallel --colsep '\t' \
        bedtools intersect \
        -a snpeff_anno/116.swapped.{1}_{2}_{3}.MODIFIER.pos.bed \
        -b /storage/zhenyingLab/huangruoshi/20211122_nh/autosome_all_neutral_25kb.apart_1kbset.bed -wa -wb \
        \| cut -f 1-3 \
        \> snpeff_anno/116.swapped.{1}_{2}_{3}.MODIFIER.NH.pos.bed \
        :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/GG_coords.txt

cat snpeff_anno/*MODIFIER.NH.pos.bed | sort -k 1,1n -k 2,2n > 116.neutral.pos.bed
```

### Deposition of coordinates of the three variant impact classes:

```
/storage/zhenyingLab/huangruoshi/genload/116.neutral.pos.bed
/storage/zhenyingLab/huangruoshi/genload/116.synonymous.pos.bed
/storage/zhenyingLab/huangruoshi/genload/116.deleterious.pos.bed
```

## Additional Note:
20240123: We have tried different neutral regions conditions, such as ensembl.cds.fa and gff.cds and repetitive region (ucsc and ensembl hardmask files). See RH_20240116.pptx
Up until this date we decided to use mRNA+ucsc, the neutral region chuck bed is in: 
```
/storage/zhenyingLab/huangruoshi/genload/116.ucsc.mRNA.neutral.auto.pos.bed
```
