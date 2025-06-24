#!/bin/bash
## Use current working directory (change working directory)
##Error file
#$ -e /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output/FST_SACE_$JOB_ID.err
## Out file
#$ -o /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output/FST_SACE_$JOB_ID.out
#$ -S /bin/bash
## Job's name
#$ -N GetFSTMatrix
#$ -l vf=12G
#$ -pe openmp 8
#$ -m e
## notification
#$ -M jabrahamavelar@gmail.com
## Modules
module load vcftools/0.1.14 
#!/bin/bash

#!/bin/bash
#$ -N map_regionB
#$ -cwd
#$ -pe openmp 8
#$ -l vf=12G
#$ -o log_regionB.out
#$ -e log_regionB.err

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

# Index reference (if not already done)
if [ ! -e "${REF}.bwt" ]; then
    bwa index "$REF"
fi
if [ ! -e "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

# Procesar cada muestra
tail -n +2 "$CSV" | while IFS=',' read -r SAMPLE _; do
    echo "Procesando muestra: $SAMPLE"

    R1="${FQDIR}/${SAMPLE}_R1_clean.fastq.gz"
    R2="${FQDIR}/${SAMPLE}_R2_clean.fastq.gz"
    UP1="${FQDIR}/${SAMPLE}_unpaired1_clean.fastq.gz"
    UP2="${FQDIR}/${SAMPLE}_unpaired2_clean.fastq.gz"
    BAM="${OUTDIR}/${SAMPLE}.bam"
    BAM_RMDUP="${OUTDIR}/${SAMPLE}.rmdup.bam"
    BAM_RG="${OUTDIR}/${SAMPLE}.rmdup.addgp.bam"
    COVERAGE_FILE="${COVDIR}/${SAMPLE}_perbase_coverage.txt"

    # Mapear
    bwa mem -M -t 8 "$REF" "$R1" "$R2" | samtools view -hbS - | samtools sort -@ 4 -o "$BAM" -
    
    # Validar BAM
    picard ValidateSamFile I="$BAM" MODE=SUMMARY O="${TMPDIR}/${SAMPLE}_bam_0_status.txt"

    # Remover duplicados
    picard MarkDuplicates INPUT="$BAM" OUTPUT="$BAM_RMDUP" METRICS_FILE="${STATDIR}/${SAMPLE}_duplicateMetrics.txt" VALIDATION_STRINGENCY=LENIENT

    # Validar post-duplicados
    picard ValidateSamFile I="$BAM_RMDUP" MODE=SUMMARY O="${TMPDIR}/${SAMPLE}_bam_1_status.txt"

    # Agregar grupos de lectura
    picard AddOrReplaceReadGroups I="$BAM_RMDUP" O="$BAM_RG" LB="$SAMPLE" PL=illumina PU="$SAMPLE" SM="$SAMPLE" VALIDATION_STRINGENCY=LENIENT

    # Validar BAM final
    picard ValidateSamFile I="$BAM_RG" MODE=SUMMARY O="${TMPDIR}/${SAMPLE}_bam_2_status.txt"

    # Indexar BAM final
    samtools index "$BAM_RG"

    # Obtener cobertura por base
    samtools depth -a "$BAM_RG" > "$COVERAGE_FILE"

    # Calcular estadísticas resumen
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
done


###HASTAQUI


# =============================
# CONFIGURACIÓN DE VARIABLES
# =============================
CSV="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE488_sinAR5.csv"
fq_dir="/mnt/Timina/lmorales/Public/ymez/data/fastq/clean"
ref_dir="/mnt/Timina/lmorales/Public/ymez/data/ref/RegionB/fasta"
ref_fasta="${ref_dir}/RegionB.fasta"
ref_short="RB"
out_dir="/mnt/Timina/lmorales/javelar/MapRegionB"
tmp="${out_dir}/tmp"
stat_dir="${out_dir}/stats"
fig_dir="${out_dir}/figs"
cov_dir="${out_dir}/coverage"
mkdir -p "$out_dir" "$tmp" "$stat_dir" "$fig_dir" "$cov_dir"

# =============================
# CARGA DE MÓDULOS
# =============================
module load bwa/0.7.4 htslib/1.2.1 gcc/5.1.0 samtools/1.9 picard/2.6.0 r/3.6.1

# =============================
# INDEXAR REFERENCIA SI ES NECESARIO
# =============================
if [ ! -f "${ref_fasta}.bwt" ]; then
    echo "Indexando referencia con bwa..."
    bwa index "$ref_fasta"
fi

if [ ! -f "${ref_fasta}.fai" ]; then
    echo "Indexando referencia con samtools..."
    samtools faidx "$ref_fasta"
fi

# =============================
# PROCESAR CADA MUESTRA
# =============================
tail -n +2 "$CSV" | cut -d',' -f1 | while read sample; do
    echo "Procesando muestra: $sample"

    # =============================
    # PASO 1: MAPEO Y ORDENAMIENTO
    # =============================
    bwa mem -M -t 10 "$ref_fasta" \
        "${fq_dir}/${sample}_R1_clean.fastq.gz" \
        "${fq_dir}/${sample}_R2_clean.fastq.gz" \
        | samtools view -hbS - \
        | samtools sort -o "${out_dir}/${sample}_${ref_short}.bam"

    # =============================
    # PASO 2: VALIDACIÓN
    # =============================
    picard ValidateSamFile \
        I="${out_dir}/${sample}_${ref_short}.bam" \
        MODE=SUMMARY \
        O="${tmp}/${sample}_${ref_short}_bam_0_status.txt"

    # =============================
    # PASO 3: MÉTRICAS DE CALIDAD
    # =============================
    picard CollectGcBiasMetrics \
        R="$ref_fasta" \
        I="${out_dir}/${sample}_${ref_short}.bam" \
        O="${stat_dir}/${sample}_${ref_short}_GCBias.txt" \
        CHART="${fig_dir}/${sample}_${ref_short}_GCBias.pdf" \
        ASSUME_SORTED=true \
        SUMMARY_OUTPUT="${stat_dir}/${sample}_${ref_short}_summary_metrics.txt" \
        VALIDATION_STRINGENCY=LENIENT

    picard MeanQualityByCycle \
        R="$ref_fasta" \
        I="${out_dir}/${sample}_${ref_short}.bam" \
        O="${stat_dir}/${sample}_${ref_short}_Qcycle.txt" \
        CHART="${fig_dir}/${sample}_${ref_short}_Qcycle.pdf" \
        VALIDATION_STRINGENCY=LENIENT

    picard QualityScoreDistribution \
        R="$ref_fasta" \
        I="${out_dir}/${sample}_${ref_short}.bam" \
        O="${stat_dir}/${sample}_${ref_short}_Qdist.txt" \
        CHART="${fig_dir}/${sample}_${ref_short}_Qdist.pdf" \
        VALIDATION_STRINGENCY=LENIENT

    # =============================
    # PASO 4: ELIMINAR DUPLICADOS
    # =============================
    picard MarkDuplicates \
        INPUT="${out_dir}/${sample}_${ref_short}.bam" \
        OUTPUT="${out_dir}/${sample}_${ref_short}.rmdup.bam" \
        METRICS_FILE="${stat_dir}/${sample}_${ref_short}_duplicateMatrix" \
        VALIDATION_STRINGENCY=LENIENT

    picard ValidateSamFile \
        I="${out_dir}/${sample}_${ref_short}.rmdup.bam" \
        MODE=SUMMARY \
        O="${tmp}/${sample}_${ref_short}_bam_1_status.txt"

    # =============================
    # PASO 5: AÑADIR GRUPOS DE LECTURA
    # =============================
    picard AddOrReplaceReadGroups \
        I="${out_dir}/${sample}_${ref_short}.rmdup.bam" \
        O="${out_dir}/${sample}_${ref_short}.rmdup.addgp.bam" \
        LB="${sample}" PL=illumina PU="${sample}" SM="${sample}" \
        VALIDATION_STRINGENCY=LENIENT

    picard ValidateSamFile \
        I="${out_dir}/${sample}_${ref_short}.rmdup.addgp.bam" \
        MODE=SUMMARY \
        O="${tmp}/${sample}_${ref_short}_bam_2_status.txt"

    samtools index "${out_dir}/${sample}_${ref_short}.rmdup.addgp.bam"

    # =============================
    # PASO 6: COBERTURA POR BASE
    # =============================
    samtools depth -aa "${out_dir}/${sample}_${ref_short}.rmdup.addgp.bam" > "${cov_dir}/${sample}_${ref_short}_depth.txt"

    # =============================
    # PASO 7: ESTADÍSTICAS DE COBERTURA
    # =============================
    awk '{print $3}' "${cov_dir}/${sample}_${ref_short}_depth.txt" | \
        Rscript -e '
            x <- scan("stdin", quiet=TRUE)
            stats <- c(mean=mean(x), sd=sd(x), median=median(x), q25=quantile(x, 0.25), q75=quantile(x, 0.75))
            write.table(t(stats), file="'${cov_dir}/${sample}_${ref_short}_coverage_stats.txt'", sep="\t", col.names=NA, quote=FALSE)
        '

    # =============================
    # LIMPIEZA
    # =============================
    if [[ -s "${out_dir}/${sample}_${ref_short}.rmdup.addgp.bam" ]]; then
        rm -f "${out_dir}/${sample}_${ref_short}.rmdup.bam"
    fi

    echo "✔️  Terminado $sample"
done
