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

vcftools --gzvcf $parentDIR/anc_inferrence/116combineVCFs/116.aO.vcf.gz \
       --keep $TXTDIR/$pop.popline.txt \
       --non-ref-ac 1 --recode \
       --out Wang2020_hglft_SYW/vcfs/$pop.NRA1.aO

bgzip Wang2020_hglft_SYW/vcfs/$pop.NRA1.aO.recode.vcf
tabix Wang2020_hglft_SYW/vcfs/$pop.NRA1.aO.recode.vcf.gz

vcftools --gzvcf Wang2020_hglft_SYW/vcfs/$pop.NRA1.aO.recode.vcf.gz \
       --freq \
       --out Wang2020_hglft_SYW/$pop.aO
```



#awk '{print$1"\t"$2-1"\t"$2"\t"$6}' WLH.aO.frq.count > WLH.aO.frq.count.bed
#awk '{print$1"\t"$2-1"\t"$2"\t"$6}' SK.aO.frq.count >SK.aO.frq.count.bed
#awk '{print$1"\t"$2-1"\t"$2"\t"$6}' YVC.aO.frq.count > YVC.aO.frq.count.bed
#awk '{print$1"\t"$2-1"\t"$2"\t"$6}' RJF.aO.frq.count > RJF.aO.frq.count.bed
#then vi to remove the header.
#colnames: chrm,start,end,derived_alleles_count

#bedtools intersect -a Wang2020_hglft_SYW/$pop.aO.frq.count.bed -b $parentDIR/wang2020_filtered/Sel_sweep_stats/hglft_genome_380ca_177970.GG6a.simple.bed -wa > Wang2020_hglft_SYW/$pop.aO.sel.count.bed
#bedtools intersect -a Wang2020_hglft_SYW/$pop.aO.frq.count.bed -b $parentDIR/wang2020_filtered/Sel_sweep_stats/hglft_genome_380ca_177970.GG6a.simple.bed -v > Wang2020_hglft_SYW/$pop.aO.nonsel.count.bed

#parallel rm $pop.{1}_{2}.count ::: sel nonsel ::: deleterious synonymous neutral
#parallel bedtools intersect -a Wang2020_hglft_SYW/$pop.aO.{1}.count.bed -b 116.aO_miyata_{2}.pos.bed \> Wang2020_hglft_SYW/$pop.{1}_{2}.count ::: sel nonsel ::: deleterious synonymous neutral

#wc -l *.*sel*count | grep -v total | awk '{print$2"\t"$1}' | sed 's/\.count//g' | sort > 4pops.region_class.num_var_sites
#parallel awk \'{sum+=\$4} END {print \"{}\\t\" sum}\' {} ::: *.*sel*count  | sed 's/\.count//g' | sort > 4pops.region_class.num_der_alleles

