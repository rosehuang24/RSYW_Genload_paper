#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --qos normal
#SBATCH -p amd-ep2
#SBATCH -J general
#SBATCH -o ./report/general.%A_%a.out
#SBATCH -e ./report/general.%A_%a.error
#SBCH --mem=10G
#SBH --array=1
#TCH --array=34-69,78-85
#SBA --array=5-28

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
ROHDIR=$parentDIR/ROH
REFDIR=$parentDIR/chicken_ref
scriptDIR=$parentDIR/useful_scripts
dataframe=$CRNTDIR/datatable/

soft=~/biosoft/
snpeff=$soft/snpEff/snpEff.jar

#pop=`head -n ${SLURM_ARRAY_TASK_ID} all.cohort.names |tail -n1 | awk '{print $1}'`
#indv=`head -n ${SLURM_ARRAY_TASK_ID}  $TXTDIR/RSYW.popline.txt |tail -n1 | awk '{print $1}'`

posDIR=$parentDIR/anc_inferrence/indv_vcfs
#this folder contains positions of aO swapped sites for all RSYW indvs, hom,het and sites
#hom: $posDIR/$indv.hom.pos.bed
#het: $posDIR/$indv.het.pos.bed
#site: $posDIR/$indv.variant_sites.pos.bed

#not run (20250806)
#count total number of alleles:
#wc -l $posDIR/*.hom.pos.bed |grep -v total | sed 's/\//\t/g;s/\./\t/g' | awk '{print$7"\t"$1}'  | sort > mid_process_file/general/total_alleles/indv.hom.count
#
#wc -l $posDIR/*.het.pos.bed |grep -v total | sed 's/\//\t/g;s/\./\t/g' | awk '{print$7"\t"$1}'  | sort > mid_process_file/general/total_alleles/indv.het.count

gerpdir=$REFDIR/GERP.relase100_58sauropsids/scores_chrms/
chickgerp=$gerpdir/chr${SLURM_ARRAY_TASK_ID}.gerp.scores.rep_removed.bed
chickgerp1=$gerpdir/chr${SLURM_ARRAY_TASK_ID}.gerp1.rep_removed.bed
chickgerp2=$gerpdir/chr${SLURM_ARRAY_TASK_ID}.gerp2.rep_remove.bed
chickNGTV=$gerpdir/chr${SLURM_ARRAY_TASK_ID}.gerp_negative.rep_removed.bed
WGgerp2=$REFDIR/GERP.relase100_58sauropsids/gerp2.WG.no_rep.bed


#=============================================
# five classes (Fig 2)
#=============================================
#parallel bedtools intersect -a $posDIR/{1}.{2}.pos.bed -b {3}.pos.bed \> mid_process_file/general/{1}-{2}-{3} :::: $TXTDIR/RSYW.popline.txt ::: hom het ::: LoF deleterious synonymous neutral gerp2 gerp226 gerp32

cd mid_process_file/general/
wc -l *-h* | grep -v total | sed 's/-/\t/g' | awk '{print$2"\t"$3"\t"$4"\t"$1}' > $dataframe/WG.indv.zyg.class.count
sed -i '1i\indv\tzyg\tclass\tcount' $dataframe/WG.indv.zyg.class.count

##
#=============================================
# GERP load
#=============================================
#parallel bedtools intersect -a $posDIR/{1}.{2}.pos.bed -b $chickgerp1 -wa -wb \| cut -f 1-3,7 \> mid_process_file/general/gerp_load/{1}-{2}-chr${SLURM_ARRAY_TASK_ID}_gerp1 :::: $TXTDIR/RSYW.popline.txt ::: hom het
#parallel bedtools intersect -a $posDIR/{1}.{2}.pos.bed -b $WGgerp2 -wa -wb \| cut -f 1-3,7 \> mid_process_file/general/gerp_load/{1}-{2}-gerp2 :::: $TXTDIR/RSYW.popline.txt ::: hom het

#cd mid_process_file/general/gerp_load

#wc -l *het-gerp2 | grep -v total | sed 's/-/\t/g' | awk '{print$2"\t"$1}' | sort > indv.het.gerp2.sites
#wc -l *hom-gerp2 | grep -v total | sed 's/-/\t/g' | awk '{print$2"\t"$1}' | sort > indv.hom.gerp2.sites
#parallel cat {1}-hom-gerp2 \|awk \'{sum+=\$4} END {print \"{1}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt | sort > indv.hom.gerp2.sum
#parallel cat {1}-het-gerp2 \|awk \'{sum+=\$4} END {print \"{1}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt | sort > indv.het.gerp2.sum


#------------------
#CDS or not
#------------------

#parallel bedtools intersect -a neutral_region-work/gff3.release110.CDS.no_rep.bed -b mid_process_file/general/gerp_load/{1}-{2}-gerp2 -wb \> mid_process_file/general/gerp_load/{1}-{2}-CDSgerp2 :::: $TXTDIR/RSYW.popline.txt ::: hom het 
#parallel bedtools intersect -a mid_process_file/general/gerp_load/{1}-{2}-gerp2  -b neutral_region-work/gff3.release110.CDS.no_rep.bed -v -wa \> mid_process_file/general/gerp_load/{1}-{2}-nonCDSgerp2 :::: $TXTDIR/RSYW.popline.txt ::: hom het 


#cd mid_process_file/general/gerp_load
#wc -l *het-CDSgerp2 | grep -v total | sed 's/-/\t/g' | awk '{print$2"\t"$1}' | sort > indv.het.CDSgerp2.sites
#wc -l *hom-CDSgerp2 | grep -v total | sed 's/-/\t/g' | awk '{print$2"\t"$1}' | sort > indv.hom.CDSgerp2.sites
#wc -l *het-nonCDSgerp2 | grep -v total | sed 's/-/\t/g' | awk '{print$2"\t"$1}' | sort > indv.het.nonCDSgerp2.sites
#wc -l *hom-nonCDSgerp2 | grep -v total | sed 's/-/\t/g' | awk '{print$2"\t"$1}' | sort > indv.hom.nonCDSgerp2.sites

#parallel cat {1}-hom-CDSgerp2 \|awk \'{sum+=\$7} END {print \"{1}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt | sort > indv.hom.CDSgerp2.sum
#parallel cat {1}-het-CDSgerp2 \|awk \'{sum+=\$7} END {print \"{1}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt | sort > indv.het.CDSgerp2.sum
#parallel cat {1}-hom-nonCDSgerp2 \|awk \'{sum+=\$4} END {print \"{1}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt | sort > indv.hom.nonCDSgerp2.sum
#parallel cat {1}-het-nonCDSgerp2 \|awk \'{sum+=\$4} END {print \"{1}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt | sort > indv.het.nonCDSgerp2.sum


#------------------
#deleterious or not
#------------------

#parallel bedtools intersect -a mid_process_file/general/gerp_load/{1}-{2}-gerp2 -b deleterious.pos.bed -wa \> mid_process_file/general/gerp_load/{1}-{2}-deleteriousgerp2 :::: $TXTDIR/RSYW.popline.txt ::: hom het
#parallel bedtools intersect -a mid_process_file/general/gerp_load/{1}-{2}-gerp2 -b LoF.pos.bed -wa \> mid_process_file/general/gerp_load/{1}-{2}-LoFgerp2 :::: $TXTDIR/RSYW.popline.txt ::: hom het

#parallel bedtools intersect -a mid_process_file/general/gerp_load/{1}-{2}-gerp2 -b synonymous.pos.bed -wa \> mid_process_file/general/gerp_load/{1}-{2}-synonymousgerp2 :::: $TXTDIR/RSYW.popline.txt ::: hom het 





#parallel bedtools intersect -a mid_process_file/general/gerp_load/{1}-{2}-gerp2 -b tmp_dir/RSYW.aO.MODERATE.pos.bed -wa \> mid_process_file/general/gerp_load/{1}-{2}-moderategerp2 :::: $TXTDIR/RSYW.popline.txt ::: hom het

#parallel cat {1}*{2}* \| awk \'{sum+=\$4} END {print \"{1}\\t{2}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt ::: hom het | sort > 1
#parallel cat {1}*{2}* \| awk \'{sum+=\$3-\$2} END {print \"{1}\\t{2}\\t\"sum}\' :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt ::: hom het | sort > 2

#paste 1 2| cut -f 1,2,3,6> indv.zyg.gerp2sum_sites
#sed -i '1i\indv\tzyg\tgerp1sum\tsites' indv.zyg.gerp2sum_sites






