#!/bin/bash

# Ruta al script base
SCRIPT="run_regionB_mapping_array.sh"

# Ruta al CSV con las muestras
CSV="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE488_sinAR5.csv"

# Contar número de muestras (sin el encabezado)
N=$(($(wc -l < "$CSV") - 1))

echo "Número de muestras: $N"

# Verificar si el script contiene el marcador para reemplazo
if ! grep -q "#\$ -t 1-" "$SCRIPT"; then
  echo "ERROR: No se encontró línea con '#$ -t 1-XXX' en $SCRIPT"
  exit 1
fi

# Crear copia temporal del script con el rango actualizado
TMP_SCRIPT="run_regionB_mapping_array_tmp.sh"
cp "$SCRIPT" "$TMP_SCRIPT"
sed -i "s/#\$ -t 1-.*/#\$ -t 1-$N/" "$TMP_SCRIPT"

# Enviar el job array a SGE
qsub "$TMP_SCRIPT"
