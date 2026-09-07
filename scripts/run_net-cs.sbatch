#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313
#SBATCH --exclusive
#SBATCH --mem-bind=local


TREE_QMC=/fs/cbcb-lab/ekmolloy/jdai123/vpg-study/software/TREE-QMC/build/tree-qmc

GT=$1
TOB=$2
OUTPREFIX=$3
TAXA=$4
NET_ID=$5
TRUE_NET=$6
MODE=$7
TIMER="/usr/bin/time"
ALPHA=$8
BETA=$9



### print all input arguments for debugging
echo "GT: $GT"
echo "TAXA: $TAXA"
echo "OUTPREFIX: $OUTPREFIX"
echo "NET_ID: $NET_ID"
echo "TRUE_NET: $TRUE_NET"
echo "MODE: $MODE"
echo "ALPHA: $ALPHA"
echo "BETA: $BETA"
echo "TOB: $TOB"

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


######### compute network ###########
NETWORK_OUTPUT="${OUTPREFIX}_alpha${ALPHA}_beta${BETA}.network.nwk"
NETWORK_USAGE_LOG="${OUTPREFIX}_alpha${ALPHA}_beta${BETA}.network.usage.log"

${TIMER} -v -o ${NETWORK_USAGE_LOG} ${TREE_QMC} --network -i ${GT} --at ${TOB} -o ${NETWORK_OUTPUT} --override
TIME_OF_NETWORK_BUILDING=$(wall_time_seconds ${NETWORK_USAGE_LOG} | awk '{print $2}')


TOTAL_TIME=$TIME_OF_NETWORK_BUILDING

CSV="${NETWORK_OUTPUT}.csv"
COMP_LINE=$(julia /fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/scripts/compare_network.jl "${NETWORK_OUTPUT}" "$TRUE_NET" false)
: > "${CSV}"
echo "taxa(round),network_id,true_network,estimated_network,method,alpha,beta,nl,e1,e2,hwcd,normalized-hwcd,wall_clock_seconds" > "${CSV}"
echo "${TAXA},${NET_ID},${TRUE_NET},${NETWORK_OUTPUT},${MODE},${ALPHA},${BETA},${COMP_LINE},${TOTAL_TIME}" >> ${CSV}

/fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/scripts/camus-dataset/run_compare_hybrid_vec.sh ${TRUE_NET} ${NETWORK_OUTPUT} ${NETWORK_OUTPUT}.hybrid_vec_comparison.csv