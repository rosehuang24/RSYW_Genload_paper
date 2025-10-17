#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p amd-ep2
#SBATCH -J neutral
#SBATCH -o ./report/neutral.%A_%a.out
#SBATCH -e ./report/neutral.%A_%a.error
#SBATCH --array=1-107
#SBATCH --mem=5G

source ~/.bash_profile
#module load /soft/modules/modulefiles/bioinfo/gatk-4.2.0.0
#module load /soft/modules/modulefiles/bioinfo/picard-2.25.1
module load vcftools/0.1.16
module load htslib/1.12

parentDIR=/storage/zhenyingLab/huangruoshi
REFDIR=/storage/zhenyingLab/huangruoshi/chicken_ref
REFfa=/storage/zhenyingLab/huangruoshi/chicken_ref/Gallus_gallus.GRCg6a.dna.toplevel.fa

CRNTDIR=$parentDIR/108
TXTDIR=/storage/zhenyingLab/huangruoshi/txt_might_be_useful
scriptDIR=$parentDIR/useful_scripts

##
#zerofour site pipeline: 
zf=$parentDIR/zerofour/zf_20230612_rerun_everything
soft=/home/zhenyingLab/huangruoshi/biosoft
indv=`head -n ${SLURM_ARRAY_TASK_ID} $TXTDIR/107.breed_indv_depth.DL.txt |tail -n1 | awk '{print $2}'`

#get coordinates of all heterozygous sites for each indv
grep -v "#" indv_vcfs/$indv.het.WG.input4.recode.vcf | awk -v OFS='\t' '{print $1,$2-1,$2}' > zerofour/indv_pos/$indv.het.pos.bed


# get whoever is in 0fold or 4fold site
bedtools intersect -a zerofour/indv_pos/$indv.het.pos.bed -b $zf/0-fold-sorted_202306.bed > zerofour/zerofold/$indv
bedtools intersect -a zerofour/indv_pos/$indv.het.pos.bed -b $zf/4-fold-sorted_202306.bed > zerofour/fourfold/$indv


