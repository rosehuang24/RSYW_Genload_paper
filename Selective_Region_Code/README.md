# Genetic load counting variants for selective region (RSYW)

Selected region for SK, YVC and WLH were detected using rank based method -- simply by multiplying rank of xpehh, fst and pi-ratio for every window. and cut-off is top 1%, so the numbers of selected windows were the same for all three breeds.

Then we extracted variant information (VCF) for RJF, YVC, SK and WLH in the seleted regions of the three breeds. In total there are twelve ```${pop}_in_${selregion}.vcf```


```
parallel bedtools intersect \
    -a ../vcfs/{1}.{2}.non-ref-ac1.recode.vcf \
    -b ../../108/RSYW_selregion/rank_CMS/{3}_rankCMS.xrf_ave.selected_regions.bed \
    -wa -wb \> rank_CMS/midfiles/{1}_in_{3}.{2}.intersect \
    ::: RJF SK YVC WLH \
    ::: low_synonymous_only.noNH modifier.2643loci moderate.gerp1 \
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
