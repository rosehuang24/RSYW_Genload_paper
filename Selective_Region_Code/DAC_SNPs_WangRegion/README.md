# Pipeline for calculating DAC per SNP

the domestic region were lifted from Wang et al. 2020, using regions overlapped by two sets of LSBL and pi ratio
This file only describe work with frequency-inferred DAC
if you want absolute count, choose ```--count2``` instead of ```--freq```

## 1. Input preparation
#for seperate population

```all.cohort.names``` contains names and number of individuals in that population/cohort.

```
RJF	9
SYW	44
SK	23
WLH	13
YVC	8
```

We need frequency file from the population vcfs, and then use the AF to calculate how many alleles should there be if any missing genotype. 
```
popfile=all.cohort.names
pop=`head -n ${SLURM_ARRAY_TASK_ID} $popfile|tail -n1 | awk '{print $1}'`
samples=`head -n ${SLURM_ARRAY_TASK_ID} $popfile|tail -n1 | awk '{print $2}'`
```

I processed this step in ```/anc_inferrence/RSYW_gGBS.sh```
I also used the anc_infer.py, the primitive one, because I need the stats swapped to determine if there is a WS or SW gene convertion.

```
SNP=116combineVCFs/116.swapped.WG.vcf.gz

vcftools --gzvcf $SNP \
        --keep $TXTDIR/$pop.popline.txt \
        --non-ref-ac 1 --max-missing 0.7 \
        --recode --out pop_vcfs/$pop.swapped.NRA1.max_missing0.7

bgzip pop_vcfs/$pop.swapped.NRA1.max_missing0.7.recode.vcf
tabix pop_vcfs/$pop.swapped.NRA1.max_missing0.7.recode.vcf.gz
```

## 2. calculate allele frequencies

```
vcftools --gzvcf $parentDIR/anc_inferrence/pop_vcfs/$pop.swapped.NRA1.max_missing0.7.recode.vcf.gz \
       --freq \
       --out  Wang2020_hglft_SYW/freq.mid.files/$pop.stat_swapped
```

## 3. process the frequency results 
Split them into selective and non selective region with three categories of impact (deleterious, synonymous and neutral)

```
#1. format transformation
awk '{print$1"\t"$2-1"\t"$2"\t"$5"\t"$6}' Wang2020_hglft_SYW/freq.mid.files/$pop.stat_swapped.frq > Wang2020_hglft_SYW/freq.mid.files/$pop.stat_swapped.frq.bed

#2. region of selection
bedtools intersect \
       -a Wang2020_hglft_SYW/freq.mid.files/$pop.stat_swapped.frq.bed \ #from last step
       -b $parentDIR/wang2020_filtered/Sel_sweep_stats/hglft_genome_380ca_177970.GG6a.simple.bed \
       -wa > Wang2020_hglft_SYW/freq.mid.files/$pop.sel.stat_swapped.frq.bed
bedtools intersect \
       -a Wang2020_hglft_SYW/freq.mid.files/$pop.stat_swapped.frq.bed \
       -b $parentDIR/wang2020_filtered/Sel_sweep_stats/hglft_genome_380ca_177970.GG6a.simple.bed \
       -v > Wang2020_hglft_SYW/freq.mid.files/$pop.nonsel.stat_swapped.frq.bed

#3. class of variant impact
parallel bedtools intersect \
       -a Wang2020_hglft_SYW/freq.mid.files/$pop.{1}.stat_swapped.frq.bed \ #from last step
       -b 116.aO_miyata_{2}.pos.bed \
       \> Wang2020_hglft_SYW/BGC/$pop.{1}_{2}.frq.bed ::: sel nonsel ::: deleterious synonymous neutral

```

## 4. use python script to calculate S, DAC and DAC_SNPs
```
parallel python3 $scriptDIR/cal_for_all.gBGC.py \
       -I Wang2020_hglft_SYW/BGC/$pop.{1}_{2}.frq.bed \ #from last step
       -O Wang2020_hglft_SYW/BGC/$pop.{1}_{2}.conversion.S.DAC.DAC_SNPs \
       -n $samples -p $pop \
       ::: sel nonsel ::: deleterious synonymous neutral
```

## 5. combine dataframe and plot with R

```
awk '{print FILENAME"\t"$0}' *S.DAC.DAC_SNPs | sed 's/\.s/\ts/g;s/\.n/\tn/g;s/\.c/\t/g;s/_/\t/g' | cut -f 1-3,7- > BGS.DOMsel.dataframe
sed -i '1i\breed\tregion\tclass\tconversion\tS\tDAC\tDAC_SNPs' BGS.DOMsel.dataframe 
```

