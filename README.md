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
      >snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.LOW.synonymous_only.bed

grep MODIFIER vcfs/116.swapped.${chrm}_${startpos}_${endpos}.SnpEff_ann.out \
      | grep -v HIGH | grep -v MODERATE | grep -v LOW | awk '{print$1"\t"$2-1"\t"$2}' \
      >  snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.MODIFIER.pos.bed

```
### 3. Further standardize the impact into deleterious, synonymous and neutral (GERP score and neutral region from ```NH_ZF_CLBL_etc/G-PHOCS_freeman_NH.md```
GERP score: https://ftp.ensembl.org/pub/release-104/compara/conservation_scores/27_sauropsids.gerp_conservation_score/README

```
parallel bedtools intersect -a snpeff_anno/116.swapped.${chrm}_${startpos}_${endpos}.MODERATE.pos.bed -b $REFDIR/GERP_score_bins/{}.bed -wa -wb \> GERP_dis/116.${chrm}_${startpos}_${endpos}.MODERATE.gerp_{} :::: GERPFILES.list
```
