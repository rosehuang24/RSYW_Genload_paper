# Code for calculating load in genome within and out of ROH region.

## 1. Extract hom/het variants or different regions

### ROH as in general

20230825 thoughts:

Although I did ROH on the swapped.vcf, I then realized that we should use the original vcf to call ROH since the segregating sites/heterozygous sites could be discarded due to variants being poor or heterozygous in GJF.
This can result in overestimation of ROH. 

Now the work is to use the ROH called from original ROH (XX.ROHs) and then extract those regions from the swapped.vcf then census for mutation load. 

```
#paths:
CRNTDIR=/storage/zhenyingLab/huangruoshi/genload
rohfile=$CRNTDIR/ROH_sheets/$indv.ROHs

vcftools --gzvcf $indvvcf --bed $rohfile \
       --recode --recode-INFO-all \
       --non-ref-ac 2 \
       --out vcfs/$indv.roh.homALT

vcftools --gzvcf $indvvcf --exclude-bed $rohfile \
        --recode --recode-INFO-all \
        --non-ref-ac 2 \
        --out vcfs/$indv.non-roh.homALT

vcftools --gzvcf $indvvcf --exclude-bed $rohfile \
        --recode --recode-INFO-all \
       --max-non-ref-ac 1 \
       --out vcfs/$indv.het
```

### Further dividing into lengths categories (short:100-300k, long:300k-1m, mega:1m+)

example rohfile: ```$CRNTDIR/ROH_sheets/mega.ROH_sheets/SK03.mega.ROHs```

example output: ```$CRNTDIR/RSYW_ROH/mega.ROHs/SK03_mega.ROH_homalt.bed```

#### Extracting ROHs categories
```
parallel --colsep '\t' awk \'\(\$2==\"{2}\"\) \&\& \(\$5\-\$4\<300001\) {print\$3\"\\t\"\$4\"\\t\"\$5}\' $CRNTDIR/ROH_sheets/short.ROH_sheets/input.4.107.ROHs \> $CRNTDIR/ROH_sheets/short.ROH_sheets/{2}.short.ROHs :::: $TXTDIR/107.breed_indv_depth.DL.txt

parallel --colsep '\t' awk \'\(\$2==\"{2}\"\) \&\& \(\$5\-\$4\>300000\) \&\& \(\$5\-\$4\<1000001\) {print\$3\"\\t\"\$4\"\\t\"\$5}\' $CRNTDIR/ROH_sheets/long.ROH_sheets/input.4.107.ROHs \> $CRNTDIR/ROH_sheets/long.ROH_sheets/{2}.long.ROHs :::: $TXTDIR/107.breed_indv_depth.DL.txt

parallel --colsep '\t' awk \'\(\$2==\"{2}\"\) \&\& \(\$5\-\$4\>1000000\) {print\$3\"\\t\"\$4\"\\t\"\$5}\' $CRNTDIR/ROH_sheets/mega.ROH_sheets/input.4.107.ROHs \> $CRNTDIR/ROH_sheets/mega.ROH_sheets/{2}.mega.ROHs :::: $TXTDIR/107.breed_indv_depth.DL.txt
```
#### Get lengths of different ROH classes

```
parallel awk \'{sum+=\$3-\$2+1} END {print \"{2}\\t\"sum}\' $CRNTDIR/ROH_sheets/{1}.ROH_sheets/{2}.{1}.ROHs \>\> RSYW_ROH/{1}.indv.length  ::: short long mega :::: $TXTDIR/RSYW.popline.txt
parallel sort RSYW_ROH/{}.indv.length \> RSYW_ROH/{}.1 ::: short long mega
parallel mv RSYW_ROH/{}.1 RSYW_ROH/{}.indv.length ::: short long mega
parallel sed -i \'1i\\indv\\t{}_length\' RSYW_ROH/{}.indv.length ::: short long mega
```

#### Extracting coordinates for variants
```
parallel bedtools intersect \
        -a vcfs/{1}.roh.homALT.recode.vcf \
        -b ROH_sheets/{2}.ROH_sheets/{1}.{2}.ROHs \
        \| awk \'{print\$1\"\\t\"\$2-1\"\\t\"\$2}\' \
        \> RSYW_ROH/{2}.ROHs/{1}_{2}.ROH_homalt.bed  \
        :::: $TXTDIR/RSYW.popline.txt ::: short long mega
```



## 2. Census the load

### for ROH in general 
```
parallel bedtools intersect \
        -a vcfs/{1}.{3}.recode.vcf \
        -b 116.{2}.pos.bed \
        \| awk \'{print\$1\"\\t\"\$2-1\"\\t\"\$2}\' \
        \> RSYW_ROH/all.ROHs/{1}.{2}.{3}.bed \
         :::: $TXTDIR/RSYW.popline.txt  ::: deleterious synonymous neutral ::: roh.homALT  non-roh.homALT het
```
### for different ROH classes

```
parallel bedtools intersect \
      -a RSYW_ROH/{2}.ROHs/{1}_{2}.ROH_homalt.bed \
      -b 116.{3}.pos.bed \
      \> RSYW_ROH/{2}.ROHs/{1}_{2}.ROH_homalt.{3} \
      :::: $TXTDIR/RSYW.popline.txt ::: short long mega ::: deleterious synonymous neutral

parallel --dryrun wc -l RSYW_ROH/{1}.ROHs/*{2} \| grep -v total \| sed \'s/ROHs\\//\\t/g\;s/_{1}/\\t/g\' \| awk \'{print\$3\"\\t\"\$1}\' \> RSYW_ROH/{1}.ROHs/{2}.{1}ROH_homalt.cat ::: short long mega ::: deleterious synonymous neutral
```

## 3. Integrate the data

```
parallel paste RSYW_ROH/{1}.ROHs/deleterious.{1}ROH_homalt.cat RSYW_ROH/{1}.ROHs/synonymous.{1}ROH_homalt.cat RSYW_ROH/{1}.ROHs/neutral.{1}ROH_homalt.cat \| cut -f 1,2,4,6 \> RSYW_ROH/{1}_homalt.indv.del.syn.neu ::: short long mega

#parallel sed -i \'1i\\indv\\tdel_rhom\\tsyn_rhom\\tneu_rhom\' {1}_homalt.indv.del.syn.neu ::: short long mega
sed -i '1i\indv\tdel_Srhom\tsyn_Srhom\tneu_Srhom' short_homalt.indv.del.syn.neu
sed -i '1i\indv\tdel_Lrhom\tsyn_Lrhom\tneu_Lrhom' long_homalt.indv.del.syn.neu
sed -i '1i\indv\tdel_Mrhom\tsyn_Mrhom\tneu_Mrhom' mega_homalt.indv.del.syn.neu

paste RSYW.delrhom_delhet_delnrhom_neurhom_neuhet_neunrhom_synrhom_synhet_synnrhom.result short_homalt.indv.del.syn.neu long_homalt.indv.del.syn.neu mega_homalt.indv.del.syn.neu  | cut -f 1-12,14-16,18-20,22-24 > RSYW_everything_withlengthsHomCount


#I forgot to plot the linear regression against the length of corresponging group, rather than total ROH lengths

##****remember!! to mannually add 0 for individuls with no mega lengths ROHS
paste RSYW.delrhom_delhet_delnrhom_neurhom_neuhet_neunrhom_synrhom_synhet_synnrhom.result short_homalt.indv.del.syn.neu long_homalt.indv.del.syn.neu mega_homalt.indv.del.syn.neu short.indv.length long.indv.length mega.indv.length | cut -f 1-12,14-16,18-20,22-24,26,28,30 > RSYW_everything_withlengthsHomCount
```
