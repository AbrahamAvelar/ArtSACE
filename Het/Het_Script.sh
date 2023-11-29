#!/bin/bash

output_file="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/OutHet2_SNP.txt"  # /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/${element}_SACE.gt.g.vcf

printf "File,Strain,Total_Variants,Het_Count,Heterozygosity\n" > "$output_file"

#vcf_file="DS002c6_SACE.gt.g.vcf"
#vcf_file=$_1

cut -d',' -f3 /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE469/SampleSheet_SACE469_V3.csv | cat | while IFS= read -r element; do
echo "Processing: $element"
vcf_file="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/${element}_SACE.gt.SNP.g.vcf"


# Calculate heterozygosity and count heterozygous sites
het_count=$(bcftools query -f '[%GT]\n' "$vcf_file" | awk '{het_count += gsub(/0\/1|1\/0/, "")} END {print het_count}')
total_variants=$(bcftools view -H "$vcf_file" | wc -l)

# Calculate heterozygosity ratio
heterozygosity=$(awk "BEGIN {print $het_count / $total_variants}")

# Store results in the output file with three columns
strain="${vcf_file%%_*}"
filename="${vcf_file##*/}"

printf "%s,%s,%d,%d,%.4f\n" "$filename" "$element" "$total_variants" "$het_count" "$heterozygosity" >> "$output_file"
done
echo "Results stored in $output_file."
