# Code for calculating load in genome within and out of ROH region.

## 1. Extract hom/het variants or different regions

### Extracting ROH as in general

20230825 thoughts:

Although I did ROH on the swapped.vcf, I then realized that we should use the original vcf to call ROH since the segregating sites/heterozygous sites could be discarded due to variants being poor or heterozygous in GJF.
This can result in overestimation of ROH. 

Now the work is to use the ROH called from original ROH (XX.ROHs) and then extract those regions from the swapped.vcf then census for mutation load. 

```
#paths:
CRNTDIR=/storage/zhenyingLab/huangruoshi/genload
TXTDIR=/storage/zhenyingLab/huangruoshi/txt_might_be_useful
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

#### Extracting coordinates of different ROH classes
```
parallel --colsep '\t' awk \'\(\$2==\"{2}\"\) \&\& \(\$5\-\$4\<300001\) {print\$3\"\\t\"\$4\"\\t\"\$5}\' $CRNTDIR/ROH_sheets/input.4.107.ROHs.LH \> $CRNTDIR/ROH_sheets/short.ROH_sheets/{2}.short.ROHs.bed :::: $TXTDIR/107.breed_indv_depth.DL.txt

parallel --colsep '\t' awk \'\(\$2==\"{2}\"\) \&\& \(\$5\-\$4\>300000\) \&\& \(\$5\-\$4\<1000001\) {print\$3\"\\t\"\$4\"\\t\"\$5}\' $CRNTDIR/ROH_sheets/input.4.107.ROHs.LH \> $CRNTDIR/ROH_sheets/long.ROH_sheets/{2}.long.ROHs.bed :::: $TXTDIR/107.breed_indv_depth.DL.txt

parallel --colsep '\t' awk \'\(\$2==\"{2}\"\) \&\& \(\$5\-\$4\>1000000\) {print\$3\"\\t\"\$4\"\\t\"\$5}\' $CRNTDIR/ROH_sheets/input.4.107.ROHs.LH \> $CRNTDIR/ROH_sheets/mega.ROH_sheets/{2}.mega.ROHs.bed :::: $TXTDIR/107.breed_indv_depth.DL.txt

```
double check if numbers of ROH match (mind the working directory)
```
wc -l *ROHs | grep -v total | grep  -v inpu | awk '{print$2"\t"$1-1}' | sort > ROH.numbers.107.noheaders
parallel wc -l {2}.ROH_sheets/{1}.{2}.ROHs :::: $TXTDIR/107.indv.DL.txt ::: short long mega > 3classes_line_count.no_headers
parallel grep {} 3classes_line_count.no_headers \| awk \'{sum+=\$1} END {print\"{}\\t\"sum}\' :::: ../../txt_might_be_useful/107.indv.DL.txt  | sort > ROH.sum.3calsses
paste ROH.numbers.107.noheaders ROH.sum.3calsses 

```

20240115
it seems that the total number of sites are different when adding vcfs of diff.lengths ROH and total ROHs, let's check if first line of bed file is omitted by vcf or bedtools. 
```
parallel sed -i \'1i\\chrm\\tstart\\tend\' $CRNTDIR/ROH_sheets/{2}.ROH_sheets/{1}.{2}.ROHs.bed :::: $TXTDIR/107.indv.DL.txt ::: short long mega

```

#### Get lengths of different ROH classes

```
parallel awk \'{sum+=\$3-\$2+1} END {print \"{2}\\t\"sum}\' $CRNTDIR/ROH_sheets/{1}.ROH_sheets/{2}.{1}.ROHs \>\> RSYW_ROH/{1}.indv.length  ::: short long mega :::: $TXTDIR/RSYW.popline.txt
parallel sort RSYW_ROH/{}.indv.length \> RSYW_ROH/{}.1 ::: short long mega
parallel mv RSYW_ROH/{}.1 RSYW_ROH/{}.indv.length ::: short long mega
parallel sed -i \'1i\\indv\\t{}_length\' RSYW_ROH/{}.indv.length ::: short long mega
```

#### Extracting variants positions in different ROH classes

```
parallel bedtools intersect \
        -a vcfs/{1}.roh.homALT.recode.vcf \
        -b ROH_sheets/{2}.ROH_sheets/{1}.{2}.ROHs.bed \
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

parallel wc -l RSYW_ROH/{1}.ROHs/*{2} \| grep -v total \| sed \'s/ROHs\\//\\t/g\;s/_{1}/\\t/g\' \| awk \'{print\$3\"\\t\"\$1}\' \> RSYW_ROH/{1}.ROHs/{2}.{1}ROH_homalt.cat ::: short long mega ::: deleterious synonymous neutral
```

## 3. Integrate the data

### for ROH in general

```
#count the number
parallel wc -l RSYW_ROH/all.ROHs/*{1}.{2}.bed \
              \| grep -v total \
              \| sed \'s/\\.{1}/\\t/g\' \
              \| awk \'{print\$2\"\\t\"\$1}\' \
              \| sort \> RSYW_ROH/all.ROHs/{1}.{2}.count \
              ::: deleterious synonymous neutral ::: roh.homALT non-roh.homALT het
#build the dataframe
cd RSYW_ROH/all.ROHs/

paste ../../RSYW.breed.indv.totalROHlength deleterious.het.count deleterious.non-roh.homALT.count deleterious.roh.homALT.count neutral.het.count neutral.non-roh.homALT.count neutral.roh.homALT.count synonymous.het.count synonymous.non-roh.homALT.count synonymous.roh.homALT.count \
       | cut -f 1-3,5,7,9,11,13,15,17,19,21 > all.ROH.cat

sed -i '1i\breed\tindv\tROHlength\tdel_het\tdel_nonROH_hom\tdel_ROH_hom\tneu_het\tneu_nonROH_hom\tneu_ROH_hom\tsyn_het\tsyn_nonROH_hom\tsyn_ROH_hom' all.ROH.cat

```

### for different ROH classes

```
parallel paste RSYW_ROH/{1}.ROHs/deleterious.{1}ROH_homalt.cat RSYW_ROH/{1}.ROHs/synonymous.{1}ROH_homalt.cat RSYW_ROH/{1}.ROHs/neutral.{1}ROH_homalt.cat \| cut -f 1,2,4,6 \> RSYW_ROH/{1}_homalt.indv.del.syn.neu ::: short long mega

#parallel sed -i \'1i\\indv\\tdel_rhom\\tsyn_rhom\\tneu_rhom\' {1}_homalt.indv.del.syn.neu ::: short long mega 
sed -i '1i\indv\tdel_Srhom\tsyn_Srhom\tneu_Srhom' short_homalt.indv.del.syn.neu
sed -i '1i\indv\tdel_Lrhom\tsyn_Lrhom\tneu_Lrhom' long_homalt.indv.del.syn.neu
sed -i '1i\indv\tdel_Mrhom\tsyn_Mrhom\tneu_Mrhom' mega_homalt.indv.del.syn.neu
```
### put those two dataframes together

```
paste all.ROHs/all.ROH.cat long_homalt.indv.del.syn.neu  mega_homalt.indv.del.syn.neu  short_homalt.indv.del.syn.neu long.indv.length  mega.indv.length  short.indv.length \
       | cut -f -12,14-16,18-20,22-24,26,28,30 \
       > RSYW.breed.indv.ROH.all.and.3classes

```
