#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --qos normal
#SBATCH -p amd-ep2
#SBATCH -J POP_DAC
#SBATCH -o ./report/POP_DAC.%A_%a.out
#SBATCH -e ./report/POP_DAC.%A_%a.error
#SBATCH --mem=20G
#SBATCH --array=1-NNNN

source ~/.bash_profile

module unload gcc
module load gcc/9.3.0
module load perl/5.34.0
module load plink/1.90
module load bcftools/1.14
module load vcftools/0.1.16
module load htslib/1.12
module load jdk/18.0.1.1 

parentDIR=/storage/zhenyingLab/huangruoshi
CRNTDIR=$parentDIR/genload_53
TXTDIR=/storage/zhenyingLab/huangruoshi/txt_might_be_useful
REFDIR=$parentDIR/chicken_ref
scriptDIR=$parentDIR/useful_scripts
dataframe=$CRNTDIR/datatable/

soft=~/biosoft/
snpeff=$soft/snpEff/snpEff.jar

indv=`head -n ${SLURM_ARRAY_TASK_ID}  $TXTDIR/POP.popline.txt |tail -n1 | awk '{print $1}'`

posDIR=$parentDIR/anc_inferrence/indv_vcfs
#this folder contains positions of aO swapped sites for all RSYW indvs, hom,het and sites
#hom: $posDIR/$indv.hom.pos.bed
#het: $posDIR/$indv.het.pos.bed
#site: $posDIR/$indv.variant_sites.pos.bed


selDIR=$parentDIR/108/RSYW_selregion/rank_CMS/
popsel=$selDIR/archive/POP_rankCMS.xrf_ave.selected_regions.merged.bed
nonsel=$parentDIR/genload_pilot/20240415_selregion/nonsel.dump.bed


#jackknife
jack=mid_process_file/art_sel/jack/
#mkdir $jack/POP
#mkdir $jack/POP/frq_files

#vcftools --gzvcf $parentDIR/anc_inferrence/pop_vcfs/POP.swapped.NRA1.max_missing0.7.recode.vcf.gz \
#       --remove-indv $indv \
#       --freq \
#       --out $jack/POP/frq_files/no_$indv


#---------------------
## DAC
#---------------------
# jackknife

## *******************************##
## Talk asfter 202412028:
# we need another gerp2 class for artificial selection
# what we discussed is that in sel region if we ecluded the 50 kb flanking region
# on both side of CDS,
# then we are removing possible signs for hitckhiking for gerp2 variants.
# therefore the gerp2 class here are merely all the sites in selected region (too few in CDS thus no care)
## *******************************##
gerpdir=$REFDIR/GERP.relase100_58sauropsids
gerp2=RSYW.WG.gerp2.pos.bed
gerp226=RSYW.WG.gerp226.pos.bed
gerp32=RSYW.WG.gerp32.pos.bed


input=$jack/POP/frq_files/no_$indv.frq

#grep -v CHR $input | awk -F ':' '($3!=0) {print $0}'  | awk -v OFS="\t" '{print $1,$2-1,$2,substr($6,3)}' > $jack/POP/frq_files/no_$indv.all.alt.freq.bed

in=$jack/POP/frq_files/no_$indv.all.alt.freq.bed

#list sites in each jackknife pooulation that sits in the selected or nonselected region

#parallel bedtools intersect -a $selDIR/archive/{1}_rankCMS.xrf_ave.selected_regions.merged.bed -b $in -wa -wb \| cut -f 4-7 \> $jack/POP/no_${indv}-in-{1}.bed ::: SK WLH YVC
#bedtools intersect -a $nonsel -b $in -wa -wb | cut -f 4-7 > $jack/POP/no_${indv}-in-nonsel.bed

# class them

parallel bedtools intersect -a  $jack/POP/no_${indv}-in-{1}.bed -b {2}.pos.bed  -wa -wb \| cut -f 1-4 \> $jack/POP/no_${indv}-in-{1}-{2}-class ::: SK WLH YVC nonsel ::: deleterious 

parallel bedtools intersect -a  $jack/POP/no_${indv}-in-{1}.bed -b  $CRNTDIR/RSYW.WG.{2}.pos.bed  -wa -wb \| cut -f 1-4 \> $jack/POP/no_${indv}-in-{1}-{2}-class ::: SK WLH YVC nonsel ::: gerp2 gerp226 gerp32




