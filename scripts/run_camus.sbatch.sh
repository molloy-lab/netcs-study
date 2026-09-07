#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313
#SBATCH --exclusive
#SBATCH --mem-bind=local

CAMUS_PATH=/fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/software/camus/camus

GT=$1

CONSTRAT_TREE=$2
OUTPREFIX=$3
TAXA=$4
NET_ID=$5
TRUE_NET=$6
METHOD=$7
RES_DIR=$8
TIMER="/usr/bin/time"

LAST_NEWICK_FILE="${OUTPREFIX}.last_extended_newick.nwk"

USAGE_LOG="${OUTPREFIX}.usage.log"

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



${TIMER} -v -o ${USAGE_LOG} ${CAMUS_PATH} -o ${OUTPREFIX} ${CONSTRAT_TREE} ${GT}
CAMUS_CSV="${OUTPREFIX}.csv"

MATCH_NEWICK_FILE="${OUTPREFIX}.matched_extended_newick.nwk"
JULIA_SELECTER=/fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/scripts/camus-dataset/camus/select_network_by_hybird_num.jl
julia ${JULIA_SELECTER} ${TRUE_NET} ${CAMUS_CSV} ${MATCH_NEWICK_FILE}

WALL_CLOCK_TIME=$(wall_time_seconds ${USAGE_LOG} | awk '{print $2}')

CSV="${RES_DIR}/CAMUS-ASTRAL-${TAXA}_network_${NET_ID}_${ILS}.level1_network.nwk.csv"
COMP_LINE=$(julia /fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/scripts/compare_network.jl "${MATCH_NEWICK_FILE}" "$TRUE_NET" false)
: > "${CSV}"
echo "taxa(round),network_id,true_network,estimated_network,method,alpha,beta,nl,e1,e2,hwcd,normalized-hwcd,wall_clock_seconds" > "${CSV}"
echo "${TAXA},${NET_ID},${TRUE_NET},${MATCH_NEWICK_FILE},${METHOD},,,${COMP_LINE},${WALL_CLOCK_TIME}" >> ${CSV}
