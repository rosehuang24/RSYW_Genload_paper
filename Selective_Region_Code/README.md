# Genetic load counting variants for selective region (RSYW)

```
parentDIR=/storage/zhenyingLab/huangruoshi
TXTDIR=/storage/zhenyingLab/huangruoshi/txt_might_be_useful
REFDIR=/storage/zhenyingLab/huangruoshi/chicken_ref
```

Selected region for SK, YVC and WLH were detected using rank based method -- simply by multiplying rank of xpehh, fst and pi-ratio for every window. and cut-off is top 1%, so the numbers of selected windows were the same for all three breeds.

Then we extracted variant information (VCF) for all individuals of RJF, YVC, SK and WLH in the seleted regions of the three breeds. In total there are 53 X 3 = 159 ``` ${indv}_in_${selregion}.vcf```

### Input preperation

```
indv=`head -n ${SLURM_ARRAY_TASK_ID}  $TXTDIR/107.breed_indv_depth.DULO.txt | tail -n1 | awk '{print $2}'`

#I keep this seperately from other steps is because these vcfs might become handy in future. 107 of them
#SBTACH -a 1-107
vcftools --gzvcf $parentDIR/anc_inferrence/116combineVCFs/116.swapped.WG.vcf.gz \
        --keep $indv
        --non-ref-ac 1 --recode \
        --out $parentDIR/anc_inferrence/indv_vcfs/$indv.swapped.WG

# to extract selregion homozygous and heterozygous
parallel vcftools --vcf $parentDIR/anc_inferrence/indv_vcfs/$indv.swapped.WG.recode.vcf \
            --bed $parentDIR/108/RSYW_selregion/rank_CMS/{1}_rankCMS.xrf_ave.selected_regions.merged.bed \
            --non-ref-ac 2 --recode \
            --out RSYW_selregion/rank_CMS/vcfs/$indv_in_{1}.hom \
            ::: SK WLH YVC

parallel vcftools --vcf $parentDIR/anc_inferrence/indv_vcfs/$indv.swapped.WG.recode.vcf \
            --bed $parentDIR/108/RSYW_selregion/rank_CMS/{1}_rankCMS.xrf_ave.selected_regions.merged.bed \
            --max-non-ref-ac 1 --recode \
            --out RSYW_selregion/rank_CMS/vcfs/$indv_in_{1}.het \
            ::: SK WLH YVC
```


### Intersect with three categories

```
parallel bedtools intersect \
        -a RSYW_selregion/rank_CMS/vcfs/$indv_in_{1}.{2}.recode.vcf \
        -b {3} \
        \| awk \'{print\$1\"\\t\"\$2-1\"\\t\"\$2}\' \
        \> RSYW_selregion/rank_CMS/midfiles/$indv_in_{1}.{2}.{4}.pos.bed \
        ::: SK WLH YVC ::: hom het \
        ::: 116.deleterious.pos.bed 116.synonymous.pos.bed 116.ucsc.cds.neutral.auto.pos.bed \
        :::+ deleterious synonymous neutral

```

Now you have each individual's polymophysms in selected regions with variant impact annotation







#============================================================================================
Then we extracted variant information (VCF) for RJF, YVC, SK and WLH in the seleted regions of the three breeds. In total there are twelve ```${pop}_in_${selregion}.vcf```


```
parallel bedtools intersect \
    -a {2} \
    -b /storage/zhenyingLab/huangruoshi/108/RSYW_selregion/rank_CMS/{3}_rankCMS.xrf_ave.selected_regions.bed \
    -wa -wb \| awk \'{print\$1\"\\t\"\$2-1\"\\t\"\$2}\' \
    \> RSYW_selregion/rank_CMS/midfiles/{1}_in_{3}.{4}.pos.bed \
    ::: RJF SK YVC WLH \
    ::: 116.deleterious.pos.bed 116.synonymous.pos.bed 116.ucsc.cds.neutral.auto.pos.bed \
    :::+ deleterious synonymous neutral \
    :::  SK YVC WLH
```





```
parallel  python3 summarize_3cate.py  \
      -I rank_CMS/midfiles/{1}_in_{3}.{2}.intersect \
      -O rank_CMS/{1}_in_{3}.{2}.3cate_var_count \
      -P /storage/zhenyingLab/huangruoshi/txt_might_be_useful/{1}.popline.txt \
      ::: RJF SK YVC WLH \
      ::: low_synonymous_only.noNH modifier.2643loci moderate.gerp1 \
      ::: SK YVC WLH
```
```
parallel paste \
      rank_CMS/{1}_in_{2}.moderate.gerp1.3cate_var_count \
      rank_CMS/{1}_in_{2}.low_synonymous_only.noNH.3cate_var_count \
      rank_CMS/{1}_in_{2}.modifier.2643loci.3cate_var_count  \| cut -f 1-3,5,6,8,9 \| grep -v het \| awk \'{print\"{2}\\t{1}\\t\"\$0}\' \
      \>  rank_CMS/{1}_in_{2}.sel.breed.indv.del.syn.neu.cat \
      ::: RJF SK YVC WLH \
      ::: SK YVC WLH
```
```
cat rank_CMS/*sel.breed.indv.del.syn.neu.cat > rank_CMS/RSYW.ranksel.breed.indv.del.syn.neu.cat
```
```
sed -i '1i\selregion\tbreed\tindv\tdel_hom\tdel_het\tsyn_hom\tsyn_het\tneu_hom\tneu_het' rank_CMS/RSYW.ranksel.breed.indv.del.syn.neu.cat 

```
