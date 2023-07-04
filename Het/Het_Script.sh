#!/bin/bash

output_file="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/OutHet.txt"
printf "File\tStrain\tTotal_Variants\tHet_Count\tHeterozygosity\n" > "$output_file"

vcf_file="DS002c6_SACE.gt.g.vcf"
vcf_file="DK002c39_SACE.gt.g.vcf"
vcf_file="XA126c1_SACE.gt.g.vcf"
vcf_file="XA124c1_SACE.gt.g.vcf"
vcf_file="XA126c5_SACE.gt.g.vcf"
vcf_file="YMX507B07_SACE.gt.g.vcf"
vcf_file="YMX005645_SACE.gt.g.vcf"
vcf_file="XA121c18_SACE.gt.g.vcf"
vcf_file="JS863c1_SACE.gt.g.vcf"
vcf_file="JS861c1_SACE.gt.g.vcf"
vcf_file="JS208c1_SACE.gt.g.vcf"
vcf_file="XA145c3_SACE.gt.g.vcf"
vcf_file="YMX005576_SACE.gt.g.vcf"


#vcf_file=$_1

# Calculate heterozygosity and count heterozygous sites
het_count=$(bcftools query -f '[%GT]\n' "$vcf_file" | awk '{het_count += gsub(/0\/1|1\/0/, "")} END {print het_count}')
total_variants=$(bcftools view -H "$vcf_file" | wc -l)

# Calculate heterozygosity ratio
heterozygosity=$(awk "BEGIN {print $het_count / $total_variants}")

# Store results in the output file with three columns
strain="${vcf_file%%_*}"
printf "%s\t%s\t%d\t%d\t%.4f\n" "$vcf_file" "$strain" "$total_variants" "$het_count" "$heterozygosity" >> "$output_file"

echo "Results stored in $output_file."
