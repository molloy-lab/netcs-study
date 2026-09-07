#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313
#SBATCH --exclusive
#SBATCH --mem-bind=local

GT=$1
TOB=$2
OUT=$3
TAXA=$4
NET_ID=$5
ILS=$6
TRUE_NET=$7
METHOD=$8

# CORES=8 ## 8 for 50 and 100 taxa but will run out memory for 200 taxa  CORES=1 for the 200taxa even 2 will run out of memory
CORES=8
ALPHA=$9
BETA=${10}
MODEL="T3"



R_SCRIPT=/fs/cbcb-lab/ekmolloy/umd-ufl-collab/tree-of-blobs-study/software/R-4.5.1-build/bin/Rscript

TINNIK_DRIVER="/fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/scripts/nanquplus_from_tob.R"
TIMER="/usr/bin/time"


OUTPUT="${OUT}.level1_network.nwk"
USAGE_LOG="${OUTPUT}.usage.log"
CSV="${OUTPUT}.csv"

wall_time_seconds() {
    for _f in "$@"; do
        # Match the "Elapsed (wall clock) time" line (ignoring extra text)
        _ts=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $NF}' "$_f")

        if [ -z "$_ts" ]; then
            echo "$_f    0.000"   # Log file didn't have the line
            continue
        fi

        IFS=: read -r _p1 _p2 _p3 <<< "$_ts"

        if [ -n "$_p3" ]; then
            _seconds=$(awk -v h="$_p1" -v m="$_p2" -v s="$_p3" \
                       'BEGIN{printf "%.3f", h*3600 + m*60 + s}')
        else
            _seconds=$(awk -v m="$_p1" -v s="$_p2" \
                       'BEGIN{printf "%.3f", m*60 + s}')
        fi

        printf '%s\t%s\n' "$_f" "$_seconds"
    done
}




# if [ ! -e "${CSV}" ]; then
${TIMER} -v -o ${USAGE_LOG} ${R_SCRIPT} --vanilla ${TINNIK_DRIVER} "$GT" "$TOB" "$OUT" "$CORES" "$ALPHA" "$BETA" "$MODEL"
# fi



WALL_CLOCK_TIME=$(wall_time_seconds ${USAGE_LOG} | awk '{print $2}')
COMP_LINE=$(julia /fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/scripts/compare_network.jl "${OUTPUT}" "$TRUE_NET" false)
: > "${CSV}"
echo "taxa(round),network_id,true_network,estimated_network,method,alpha,beta,nl,e1,e2,hwcd,normalized-hwcd,wall_clock_seconds" > "${CSV}"
echo "${TAXA},${NET_ID},${ILS},${TRUE_NET},${OUTPUT},${METHOD},${ALPHA},${BETA},${COMP_LINE},${WALL_CLOCK_TIME}" >> ${CSV}
