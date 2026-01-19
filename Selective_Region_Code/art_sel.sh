#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --qos normal
#SBATCH -p amd-ep2
#SBATCH -J art_sel
#SBATCH -o ./report/art_sel.%A_%a.out
#SBATCH -e ./report/art_sel.%A_%a.error
#SBAT --mem=10G
#SBATCH --array=1-4
#STCH --array=1-9

#pop=RJF
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
#pop=WLH

posDIR=$parentDIR/anc_inferrence/indv_vcfs
#this folder contains positions of aO swapped sites for all RSYW indvs, hom,het and sites
#hom: $posDIR/$indv.hom.pos.bed
#het: $posDIR/$indv.het.pos.bed
#site: $posDIR/$indv.variant_sites.pos.bed


#selDIR=$parentDIR/108/RSYW_selregion/rank_CMS/
#popsel=$selDIR/archive/${pop}_rankCMS.xrf_ave.selected_regions.merged.bed
#nonsel=$parentDIR/genload_pilot/20240415_selregion/nonsel.dump.bed

## 8*******************************##
## Talk asfter 202412028:
# we need another gerp2 class for artificial selection
# what we discussed is that in sel region if we ecluded the 50 kb flanking region 
# on both side of CDS, 
# then we are removing possible signs for hitckhiking for gerp2 variants.
# therefore the gerp2 class here are merely all the sites in selected region (too few in CDS thus no care)
## 8*******************************##

selcds=mid_process_file/art_sel/denominators/$pop.cds.bed
nonselcds=mid_process_file/art_sel/denominators/nonsel.cds.bed
selnoncds=mid_process_file/art_sel/denominators/$pop.noncds.bed
nonselnoncds=mid_process_file/art_sel/denominators/nonsel.noncds.bed


# ==========================================
# DAC_SNPs
# ==========================================

#---------------------------
## jackknife preperation
#---------------------
jack=$CRNTDIR/mid_process_file/art_sel/jack/private_sel
#first run DAC_jackknife.sh
#**
# very IMPORTANT: check WLH.sh, the popline file should be the one with 13 indvs.
#**
#parallel sed \'s/POP/{1}/g\;s/NNNN/{2}/g\' DAC_jackknife.sh \> {1}.sh ::: RJF YVC SK WLH :::+ 9 8 23 13
#parallel sbatch {1}.sh ::: RJF YVC SK WLH



#---------------------
# process data
#---------------------
#remember, I did not count th DAC here, please multiply the total_freq by number of chromosomes when plotting.


POP=`head -n ${SLURM_ARRAY_TASK_ID} all.cohort.names |tail -n1 | awk '{print $1}'`
cd $jack/$POP
parallel -q awk '{sum+=$4} END {OFMT="%.0f"; print FILENAME, sum}' {} ::: *-class | sed 's/-/\t/g' | awk -v OFS='\t' -v i=$POP '{print i,$1,$3,$4,$6}' > $jack/$POP.totalfreq

#---------------------
# data tidy, no array
#---------------------
#cd $jack
#cat *freq > $dataframe/jack.breed.indv.region.class.totalfreq.private_sel
#sed -i '1i\breed\tjack_indv\tregion\tclass\ttotalfreq'  $dataframe/jack.breed.indv.region.class.totalfreq.private_sel

#wc -l */*-class | sed 's/-/\t/g;s/\//\t/g' | awk -v OFS='\t' '{print$2,$3,$5,$6,$1}'  | grep -v total > $dataframe/jack.breed.indv.region.class.site.private_sel
#sed -i '1i\breed\tjack_indv\tregion\tclass\tsite' $dataframe/jack.breed.indv.region.class.site.private_sel






