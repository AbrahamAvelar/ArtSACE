#!/bin/bash
#$ -e /mnt/Timina/lmorales/javelar/MapRegionB/MapRB_$TASK_ID.err
#$ -o /mnt/Timina/lmorales/javelar/MapRegionB/MapRB_$TASK_ID.out
#$ -S /bin/bash
#$ -N MapRegionB
#$ -l vf=12G
#$ -pe openmp 8
#$ -m e
#$ -M jabrahamavelar@gmail.com
#$ -cwd
#$ -t 1-488

# Carga de módulos
module load bwa/0.7.4 htslib/1.2.1 gcc/5.1.0 samtools/1.9 picard/2.6.0 r/3.6.1

# Variables
REF="/mnt/Timina/lmorales/Public/ymez/data/ref/RegionB/fasta/RegionB.fasta"
CSV="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE488_sinAR5.csv"
FQDIR="/mnt/Timina/lmorales/Public/ymez/data/fastq/clean"
OUTDIR="/mnt/Timina/lmorales/javelar/MapRegionB"
TMPDIR="${OUTDIR}/tmp"
STATDIR="${OUTDIR}/stats"
COVDIR="${OUTDIR}/coverage"

mkdir -p "$OUTDIR" "$TMPDIR" "$STATDIR" "$COVDIR"

# Indexar referencia si no existe
if [ ! -e "${REF}.bwt" ]; then
    bwa index "$REF"
fi
if [ ! -e "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

# Obtener muestra correspondiente al task ID
SAMPLE=$(tail -n +2 "$CSV" | sed -n "${SGE_TASK_ID}p" | cut -d',' -f1)

echo "Procesando muestra: $SAMPLE (Task ID: $SGE_TASK_ID)"

R1="${FQDIR}/${SAMPLE}_R1_clean.fastq.gz"
R2="${FQDIR}/${SAMPLE}_R2_clean.fastq.gz"
UP1="${FQDIR}/${SAMPLE}_unpaired1_clean.fastq.gz"
UP2="${FQDIR}/${SAMPLE}_unpaired2_clean.fastq.gz"
BAM="${OUTDIR}/${SAMPLE}.bam"
BAM_RMDUP="${OUTDIR}/${SAMPLE}.rmdup.bam"
BAM_RG="${OUTDIR}/${SAMPLE}.rmdup.addgp.bam"
COVERAGE_FILE="${COVDIR}/${SAMPLE}_perbase_coverage.txt"

# Mapeo
bwa mem -M -t 8 "$REF" "$R1" "$R2" | samtools view -hbS - | samtools sort -@ 4 -o "$BAM" -

# Validaciones y QC
picard ValidateSamFile I="$BAM" MODE=SUMMARY O="${TMPDIR}/${SAMPLE}_bam_0_status.txt"
picard MarkDuplicates INPUT="$BAM" OUTPUT="$BAM_RMDUP" METRICS_FILE="${STATDIR}/${SAMPLE}_duplicateMetrics.txt" VALIDATION_STRINGENCY=LENIENT
picard ValidateSamFile I="$BAM_RMDUP" MODE=SUMMARY O="${TMPDIR}/${SAMPLE}_bam_1_status.txt"
picard AddOrReplaceReadGroups I="$BAM_RMDUP" O="$BAM_RG" LB="$SAMPLE" PL=illumina PU="$SAMPLE" SM="$SAMPLE" VALIDATION_STRINGENCY=LENIENT
picard ValidateSamFile I="$BAM_RG" MODE=SUMMARY O="${TMPDIR}/${SAMPLE}_bam_2_status.txt"
samtools index "$BAM_RG"

# Cobertura por base
samtools depth -a "$BAM_RG" > "$COVERAGE_FILE"

# Estadísticas resumen
Rscript -e "
cov <- read.table('$COVERAGE_FILE')[,3];
stats <- data.frame(
  sample = '$SAMPLE',
  mean = mean(cov),
  sd = sd(cov),
  median = median(cov),
  IQR1 = quantile(cov, 0.25),
  IQR3 = quantile(cov, 0.75)
);
write.table(stats, file='${STATDIR}/summary_coverage_stats.tsv', sep='\t', col.names=!file.exists('${STATDIR}/summary_coverage_stats.tsv'), row.names=FALSE, append=TRUE)
"

# Limpieza
if [[ -s "$BAM_RG" ]]; then
    rm -f "$BAM" "$BAM_RMDUP"
fi
