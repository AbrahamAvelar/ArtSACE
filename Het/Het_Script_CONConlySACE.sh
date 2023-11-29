#!/bin/bash
module load bcftools/1.9

output_file="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Output_heterozygosity_CONC.onlychr_SACE.g.vcf.gz.txt"
printf "File,Strain,Total_Variants,Het_Count,Heterozygosity\n" > "$output_file"

cut -d',' -f1 /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/SampleSheet_SACE467_V8.csv | head | cat | while IFS= read -r element; do

echo "Processing: $element"
cp "/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/${element}_CONC.gt.SNP_onlychr_SACE.g.vcf.gz" .
gunzip "${element}_CONC.gt.SNP_onlychr_SACE.g.vcf.gz"
ls -alth | head
vcf_file="${element}_CONC.gt.SNP_onlychr_SACE.g.vcf"

# Calculate heterozygosity and count heterozygous sites
het_count=$(bcftools query -f '[%GT]\n' "$vcf_file" | awk '{het_count += gsub(/0\/1|1\/0/, "")} END {print het_count}')
total_variants=$(bcftools view -H "$vcf_file" | wc -l)

# Calculate heterozygosity ratio
heterozygosity=$(awk "BEGIN {print $het_count / $total_variants}")

# Store results in the output file with three columns
strain="${vcf_file%%_*}"
filename="${vcf_file##*/}"

printf "%s,%s,%d,%d,%.4f\n" "$filename" "$element" "$total_variants" "$het_count" "$heterozygosity" >> "$output_file"
rm "${element}_CONC.gt.SNP_onlychr_SACE.g.vcf"
done
echo "Results stored in $output_file."
