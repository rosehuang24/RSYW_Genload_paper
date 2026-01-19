#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH -J ROH
#SBATCH -o ./report/ROH.out.%A_%a
#SBATCH -e ./report/ROH.error.%A_%a
#SBATCH --array=1-53
#SBT --mem=8G
source ~/.bash_profile

#module load gcc/10.2.0
module load vcftools/0.1.16
module load plink/1.90

parentDIR=/storage/zhenyingLab/huangruoshi
CRNTDIR=$parentDIR/genload_53
statDIR=$CRNTDIR/forstats
TXTDIR=/storage/zhenyingLab/huangruoshi/txt_might_be_useful
REFDIR=$parentDIR/chicken_ref
dataframe=$CRNTDIR/datatable/
soft=/soft/modules/modulefiles/bioinfo

annovar=/home/zhenyingLab/huangruoshi/biosoft/annovar
popdecay=/home/zhenyingLab/huangruoshi/biosoft/PopLDdecay-3.41/bin
ibd=/home/zhenyingLab/huangruoshi/biosoft/ibdseq.r1206.jar


#SNPVCF=$CRNTDIR/snp.input.3.py0.9missing.vcf
#SNPVCF=snp.input.4.masked_removed.recode.vcf.gz
#INDELVCF=$CRNTDIR/input.3.INDEL.header_vt0.9_missing.vcf
#SNPVCF=$parentDIR/108/RSYW.input4.recode.vcf
 
## ==========================
# BCFtools/RoH
## ==========================
mid=$CRNTDIR/mid_process_file/ROH/

#pop=`head -n ${SLURM_ARRAY_TASK_ID} all.cohort.names |tail -n1 | awk '{print $1}'`

#bcftools +fill-tags $CRNTDIR/pop_vcfs/$pop.input.4.NRA1.recode.vcf.gz -o $CRNTDIR/pop_vcfs/$pop.input.4.NRA1.AF_updated.vcf -- -t AC,AN,AF

#prepare frequency files
#bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' $CRNTDIR/pop_vcfs/$pop.input.4.NRA1.AF_updated.vcf | bgzip -c > $mid/$pop.freqs.tab.gz
#tabix -s1 -b2 -e2 $mid/$pop.freqs.tab.gz
#bcftools roh $CRNTDIR/pop_vcfs/$pop.input.4.NRA1.AF_updated.vcf --AF-file $mid/$pop.freqs.tab.gz -o $mid/$pop.bcftools.roh.default.params 

#no array, cmd or slurm both fine
#detial process
#cd $mid
#cat YVC.bcftools.roh.default.params RJF.bcftools.roh.default.params SK.bcftools.roh.default.params WLH.bcftools.roh.default.params  | awk '($6>499999)' > RSYW.bcftools.roh500k
#parallel awk \'\(\$2==\"{}\"\) {print\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}\'  RSYW.bcftools.roh500k \> {}.ROH.bed :::: $TXTDIR/RSYW.popline.txt



# ------------------#
# to find cds proportion in each segment
# ------------------#
CDS=$CRNTDIR/neutral_region-work/gff3.release110.CDS.no_rep.bed
REP=$REFDIR/repetitive_auto.bed
CDSflanking=$CRNTDIR/neutral_region-work/gff3.release110.CDS.50kbflanking.inrange.bed
indv=`head -n ${SLURM_ARRAY_TASK_ID}  $TXTDIR/RSYW.popline.txt |tail -n1 | awk '{print $1}'`

mkdir $mid/single_segments/$indv
cd $mid/single_segments/$indv

parallel --colsep '\t' echo -e \"{1}\\t{2}\\t{3}\" \> {1}-{2}-{3}.raw.bed :::: $mid/$indv.ROH.bed

parallel --colsep '\t' bedtools subtract -a {1}-{2}-{3}.raw.bed -b $REP \> {1}-{2}-{3}.ROH.no_rep.bed :::: $mid/$indv.ROH.bed

parallel --colsep '\t' rm {1}-{2}-{3}.raw.bed :::: $mid/$indv.ROH.bed

parallel --colsep '\t' bedtools intersect -a {1}-{2}-{3}.ROH.no_rep.bed -b $CDS \> {1}-{2}-{3}.CDS.no_rep.bed :::: $mid/$indv.ROH.bed

parallel -q awk '{sum+=$3-$2} END {print FILENAME"\t"sum}'  {} ::: *ROH.no_rep.bed | sed 's/.ROH.no_rep.bed//g' | sort | awk ' { if ($2 == "") { $2 = 0 } print$0}' | sed 's/ /\t/g' > ../tmpdata/$indv.ROH.no_rep.length
parallel -q awk '{sum+=$3-$2} END {print FILENAME"\t"sum}'  {} ::: *CDS.no_rep.bed | sed 's/.CDS.no_rep.bed//g' | sort | awk ' { if ($2 == "") { $2 = 0 } print$0}' | sed 's/ /\t/g' > ../tmpdata/$indv.CDS.no_rep.length

cd ../tmpdata/
paste $indv.CDS.no_rep.length $indv.ROH.no_rep.length | awk '{ if ($1!=$3) {print$0"\t-"} print$1"\t"$2"\t"$4}' > $indv.CDS.ROH.no_rep

#cmd, no slurm
#parallel -q awk '{print FILENAME"\t"$0}' {} ::: *.CDS.ROH.no_rep | sed 's/.CDS.ROH.no_rep//g;s/-/\t/g' | awk -v OFS="\t" '{print$1,$2,$3,$4,$4-$3,$6,$5}' | awk '{print$1"\t"$2"-"$3"-"$4"\t"$5"\t"$6"\t"$7}' > ../indv.segment.repROH.norepROH.CDS
#sed -i '1i\indv\tsegment\trepROH\tnorepROH\tCDS' ../indv.segment.repROH.norepROH.CDS


#parallel -q awk '{print FILENAME"\t"$0}' {} ::: *.CDS.ROH.no_rep | sed 's/.CDS.ROH.no_rep//g;s/-/\t/g' | awk -v OFS="\t" '($2<6){print$1,$2,$3,$4,$4-$3,$6,$5}' | awk '{print$1"\t"$2"\t"$2"-"$3"-"$4"\t"$5"\t"$6"\t"$7}' > ../indv.chrm1-5.segment.repROH.norepROH.CDS
#sed -i '1i\indv\tchrm\tsegment\trepROH\tnorepROH\tCDS' ../indv.chrm1-5.segment.repROH.norepROH.CDS

## ==========================
# plot for paper
## ==========================
# 20250730: let's try bcftools 500kb, 1Mb, and 2Mb+ as three categories, short-medium-long
work=$mid/analysis
cate=$mid/categories


#---
#categorying
#indv=`head -n ${SLURM_ARRAY_TASK_ID}  $TXTDIR/RSYW.popline.txt |tail -n1 | awk '{print $1}'`
#
#awk '($4>499999) && ($4<1000000)' $mid/$indv.ROH.bed | cut -f 1-3 > $cate/$indv.short.ROH.with_rep.bed
#awk '($4>999999) && ($4<2000000)' $mid/$indv.ROH.bed | cut -f 1-3 > $cate/$indv.medium.ROH.with_rep.bed
#awk '($4>1999999)' $mid/$indv.ROH.bed | cut -f 1-3 > $cate/$indv.long.ROH.with_rep.bed

#parallel bedtools subtract -a $cate/$indv.{}.ROH.with_rep.bed -b $REP \| sort -k 1,1n -k 2,2n \|  bedtools merge -i - \> $cate/$indv.{}.no_rep.bed ::: short long medium 
#bedtools subtract -a $REFDIR/auto.no_rep.bed -b $mid/$indv.ROH.bed | sort -k 1,1n -k 2,2n |  bedtools merge -i - > $cate/$indv.nonROH.no_rep.bed


#----
#basic info: Fig.3a

#no slurm, cmd
#cd $cate 
#parallel -q awk '{sum+=$3-$2} END {print FILENAME"\t"sum}' {} ::: *.no_rep.bed | sed 's/.no_rep.bed//g;s/\./\t/g' | awk ' { if ($3 == "") { $3 = 0 } print$0}' | sed 's/ /\t/g' > ROH.indv.category.no_rep_length
#sed -i '1i\indv\tcategory\tlength' ROH.indv.category.no_rep_length
#cp ROH.indv.category.no_rep_length $dataframe


#----
#CDS proportion for  Fig.3b

## *****
## * do not start unless you finished the previous data tidy cmd
## *****
#cd $cate
#parallel bedtools intersect -a  $indv.{}.no_rep.bed -b $CDS \| sort -k 1,1n -k 2,2n \| bedtools merge -i - \> $indv.{}.CDS.no_rep.bed ::: short long medium nonROH
#parallel bedtools subtract -a  $indv.{}.no_rep.bed -b $CDSflanking \| sort -k 1,1n -k 2,2n \| bedtools merge -i - \>  $indv.{}.nonCDS.no_rep.bed ::: short long medium nonROH

#--no slurm,
#parallel awk \'{sum+=\$3-\$2} END {print \"{1}\\t{2}\\t\"sum}\' $cate/{1}.{2}.CDS.no_rep.bed :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW50.popline.txt ::: short medium long nonROH | awk '{ if ($3 =="") {$3=0} print$0}' | sed 's/ /\t/g' >  $dataframe/ROH.indv.cate.CDS.no_rep_length
#sed -i '1i\indv\tcategory\tCDS'  $dataframe/ROH.indv.cate.CDS.no_rep_length

#--------
# for chrm 1-5  SuppFig
# ---------
#cd $work/chrm1-5
#parallel  awk \'\(\$1==\"{2}\"\)\' $cate/$indv.{1}.no_rep.bed \> $indv.{1}.chrm{2}.cate.no_rep.bed ::: short long medium nonROH ::: 1 2 3 4 5
#
#parallel  awk \'\(\$1==\"{2}\"\)\' $cate/$indv.{1}.CDS.no_rep.bed \> $indv.{1}.chrm{2}.CDS.no_rep.bed ::: short long medium nonROH ::: 1 2 3 4 5

#--no slurm, cmd
#ROH cates
#parallel awk \'{sum+=\$3-\$2} END {print \"{1}\\t{2}\\t{3}\\t\"sum}\' {1}.{2}.chrm{3}.cate.no_rep.bed  :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt ::: short medium long nonROH ::: 1 2 3 4 5 | awk '{ if ($4 =="") {$4=0} print$0}' | sed 's/ /\t/g' > $dataframe/ROH.indv.category.cate.chrm1-5.norep.length
#sed -i '1i\indv\tcategory\tchrm\tlength' $dataframe/ROH.indv.category.cate.chrm1-5.norep.length

#cds
#parallel awk \'{sum+=\$3-\$2} END {print \"{1}\\t{2}\\t{3}\\t\"sum}\' {1}.{2}.chrm{3}.CDS.no_rep.bed  :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.popline.txt ::: short medium long nonROH ::: 1 2 3 4 5 | awk '{ if ($4 =="") {$4=0} print$0}' | sed 's/ /\t/g' > $dataframe/ROH.indv.category.CDS.chrm1-5.norep.length
#sed -i '1i\indv\tcategory\tchrm\tCDS' $dataframe/ROH.indv.category.CDS.chrm1-5.norep.length

#----
# HDR,: Fig.3cd
posDIR=$parentDIR/anc_inferrence/indv_vcfs
#hom=$posDIR/$indv.hom.pos.bed
#homClass=$CRNTDIR/mid_process_file/general/$indv-hom-<class>
#cd $CRNTDIR
#no. of hom vars
#parallel bedtools intersect -a $CRNTDIR/mid_process_file/general/$indv-hom-{1} -b $cate/$indv.{2}.no_rep.bed  \> $work/HDR/$indv-{1}-{2}-homcount ::: deleterious neutral synonymous LoF gerp2 gerp226 gerp32 ::: short long medium nonROH


#no slurm
#cd $work/HDR/
#wc -l *homcount | grep -v total | sed 's/-/ /g' | awk -v OFS='\t' '{print $2,$3,$4,$1}' > $dataframe/ROH.indv.category.class.homcount
#sed -i '1i\indv\tclass\tcategory\tcount' $dataframe/ROH.indv.category.class.homcount

# length for calculating density
# no slurm
#length files already at:$cate
#cate/$indv.{}.CDS.no_rep.bed ::: short long medium nonROH
#cate/$indv.{}.nonCDS.no_rep.bed ::: short long medium nonROH

#-CDS
# the same comd for Fig 3b WG CDS proportion
#parallel awk \'{sum+=\$3-\$2} END {print \"{1}\\t{2}\\t\"sum}\' $cate/{1}.{2}.CDS.no_rep.bed :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW50.popline.txt ::: short medium long nonROH | awk '{ if ($3 =="") {$3=0} print$0}' | sed 's/ /\t/g' >  $dataframe/ROH.indv.cate.CDS.no_rep_length
#sed -i '1i\indv\tcategory\tCDS'  $dataframe/ROH.indv.cate.CDS.no_rep_length

#-nonCDS with flanking
#parallel awk \'{sum+=\$3-\$2} END {print \"{1}\\t{2}\\t\"sum}\' $cate/{1}.{2}.nonCDS.no_rep.bed :::: /storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW50.popline.txt ::: short medium long nonROH | awk '{ if ($3 =="") {$3=0} print$0}' | sed 's/ /\t/g' >  $dataframe/ROH.indv.cate.nonCDS.no_rep_length
#sed -i '1i\indv\tcategory\tnonCDS'  $dataframe/ROH.indv.cate.nonCDS.no_rep_length







