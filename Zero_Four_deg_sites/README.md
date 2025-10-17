# 1. Preparation
Generate CDS seq

```
grep "^>" $REFDIR/Gallus_gallus.GRCg6a.cds.all.fa | awk -F ":" '($3<29) {print$3"\t"$4-1"\t"$5}' | sort -k1,1 -k2,2n| bedtools merge > CDS.auto.from.ensembl.bed

bedtools getfasta -fi $REFDIR/Gallus_gallus.GRCg6a.dna.toplevel.fa -bed CDS.auto.from.ensembl.bed -fo CDS.from.ensembl.ffn

```

# 2. Create substitution for every site
...by making a fake vcf
```
python3 fakemuta.py CDS.from.ensembl.ffn fake_mutation.vcf
sort -V fake_mutation.vcf | uniq > fakemuta_uniq.vcf
sed -i '1i\##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' fakemuta_uniq.vcf
```

# 3. Annotate the effects
Find the ones with three missense or synonymous

```
java -Xmx40g -jar $snpeff -c $soft/snpEff/snpEff.config -v GRCg6a.99 fakemuta_uniq.vcf > output_ann.vcf
python3 fold_categorizing.py output_ann.vcf four-fold.coords zero-fold.coords all_others-folds.coords


sort four-fold.coords > 4-fold-sorted_202306.coords
sort zero-fold.coords > 0-fold-sorted_202306.coords
```
# 4. Downstream analysis
Details in ```job.sh```. Work related with neutral sites can be found in ..
