#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --qos normal
#SBATCH -p amd-ep2
#SBATCH -J rxy
#SBATCH -o ./report/Rxy.%A_%a.out
#SBATCH -e ./report/Rxy.%A_%a.error
#SBATCH --mem=20G
#SBH --array=1
#TCH --array=34-69,78-85
#SBACH --array=1-4

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

#mind WLH, 10 or 13
#pop=`head -n ${SLURM_ARRAY_TASK_ID} all.cohort.names |tail -n1 | awk '{print $1}'`
#indv=`head -n ${SLURM_ARRAY_TASK_ID}  $TXTDIR/RSYW.popline.txt |tail -n1 | awk '{print $1}'`
pop=WLH13
posDIR=$parentDIR/anc_inferrence/indv_vcfs
#this folder contains positions of aO swapped sites for all RSYW indvs, hom,het and sites
#hom: $posDIR/$indv.hom.pos.bed
#het: $posDIR/$indv.het.pos.bed
#site: $posDIR/$indv.variant_sites.pos.bed

VCF=/storage/zhenyingLab/huangruoshi/anc_inferrence/RSYW.aO.swapped.recode.vcf
#
vcftools --vcf $VCF \
	--keep $TXTDIR/$pop.popline.txt \
	--freq2 \
	--out  mid_process_file/rxy/$pop.aO.swapped

#

cd mid_process_file/rxy/
mv $pop.aO.swapped.frq WLH.aO.swapped.frq
paste RJF.aO.swapped.frq YVC.aO.swapped.frq SK.aO.swapped.frq WLH.aO.swapped.frq | grep -v CHR | awk -v OFS='\t' '{print $1,$2-1,$2,$6,$12,$18,$24}' > site.ALT.freq.RSYW.bed

parallel bedtools intersect -a site.ALT.freq.RSYW.bed -b $CRNTDIR/{}.pos.bed -wa \> $dataframe/RXY/ALT.freq.RSYW.{}.bed ::: deleterious neutral synonymous gerp2 gerp226 gerp32 LoF

#### ---------
# IMPORTANT#
# the GERP 2 in artificially selected region should contain both coding and non-coding (atleast it is so in DAC_SNP analysis)
# it is possible that the RXY gerp2 was non-coding only, I could not find the cmd....



##To do jackknifing
#See Rsy.R 
#Apparently you can jackknife in R

##old:
#for each population, indv-jackked frequency file sits in mid_process_file/art_sel/jack/freq_file/
#cmd can be found in art_sel.sh, where I sed NNNN with population names to jack away individuals in each population. 






