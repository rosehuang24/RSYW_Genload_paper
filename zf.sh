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
angDIR=/soft/bio/angsd-0.935/angsd
CRNTDIR=$parentDIR/108
TXTDIR=/storage/zhenyingLab/huangruoshi/txt_might_be_useful
scriptDIR=$parentDIR/useful_scripts
OUTDIR=$CRNTDIR/wtheta/

#zerofour site pipeline:
zf=$parentDIR/zerofour/zf_20230612_rerun_everything

#soft=/soft/modules/modulefiles/bioinfo
soft=/home/zhenyingLab/huangruoshi/biosoft

indv=`head -n ${SLURM_ARRAY_TASK_ID} $TXTDIR/107.breed_indv_depth.DL.txt |tail -n1 | awk '{print $2}'`

grep -v "#" indv_vcfs/$indv.het.WG.input4.recode.vcf | awk -v OFS='\t' '{print $1,$2-1,$2}' > zerofour/indv_pos/$indv.het.pos.bed
#
bedtools intersect -a zerofour/indv_pos/$indv.het.pos.bed -b $zf/0-fold-sorted_202306.bed > zerofour/zerofold/$indv
bedtools intersect -a zerofour/indv_pos/$indv.het.pos.bed -b $zf/4-fold-sorted_202306.bed > zerofour/fourfold/$indv
#use cmd to generate datatable.

#Neutral Het
mkdir zerofour/tmp

cool=$parentDIR/genload/neutral_region-work/non_CDS_flanking_no_rep.bed

#split chrms
awk -v chrm=${SLURM_ARRAY_TASK_ID} '($1==chrm)' $cool > zerofour/tmp/chr${SLURM_ARRAY_TASK_ID}.cool1.bed

cool1=zerofour/tmp/chr${SLURM_ARRAY_TASK_ID}.cool1.bed

cool2=$REFDIR/GERP.relase100_58sauropsids/scores_chrms/chr${SLURM_ARRAY_TASK_ID}.gerp_negative.rep_removed.bed

bedtools intersect -a $cool1 -b $cool2 | cut -f 1-3 | sort -k 1,1n -k 2,2n | bedtools merge -i - > zerofour/tmp/chr${SLURM_ARRAY_TASK_ID}.neutral.bed

parallel awk \'\(\$1=={}\)\' zerofour/indv_pos/$indv.het.pos.bed \> zerofour/tmp/$indv.chr{}.bed ::: {1..28}

parallel bedtools intersect -a zerofour/tmp/$indv.chr{}.bed -b zerofour/tmp/chr{}.neutral.bed \> zerofour/tmp/$indv.chr{}.het.neutral.pos.bed ::: {1..28}

cd zerofour/tmp/
cat $indv.*.het.neutral.pos.bed | sort -k 1,1n -k 2,2n > ../neutral/$indv.neutral.het
