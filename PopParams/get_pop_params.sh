module load vcftools/0.1.14
for popfile in *.txt; # puede ser también csv; 
do 
echo $popfile;     
population="${popfile%.txt}"; 
echo $population;  

vcftools --vcf ../../Admixture/Matrix_SNPs_SACE_from_SACE_gt_onlySNPs_filtered_missing_10_biallelic  --keep "$popfile" --site-pi --out "${population}_results";

vcftools --vcf ../../Admixture/Matrix_SNPs_SACE_from_SACE_gt_onlySNPs_filtered_missing_10_biallelic  --keep "$popfile" --window-pi 10000 --out "${population}_Window_results";

vcftools --vcf ../../Admixture/Matrix_SNPs_SACE_from_SACE_gt_onlySNPs_filtered_missing_10_biallelic  --keep "$popfile"  --TajimaD 10000 --out "${population}_TDw_results";

vcftools --vcf ../../Admixture/Matrix_SNPs_SACE_from_SACE_gt.vcf  --keep "$popfile" --site-pi --out "${population}_results_AllMatrix";

vcftools --vcf ../../Admixture/Matrix_SNPs_SACE_from_SACE_gt.vcf  --keep "$popfile" --window-pi 10000 --out "${population}_Window_results_AllMatrix";

vcftools --vcf ../../Admixture/Matrix_SNPs_SACE_from_SACE_gt.vcf  --keep "$popfile"  --TajimaD 10000 --out "${population}_TDw_results_AllMatrix";

done
