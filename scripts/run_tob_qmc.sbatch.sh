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


TREE_QMC=/fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/software/TREE-QMC/build/tree-qmc

GT=$1
ASTRAL_TREE=$2
OUTPREFIX=$3
TAXA=$4
NET_ID=$5
TRUE_NET=$6
MODE=$7
TIMER="/usr/bin/time"
ALPHA=$8
BETA=$9
TRUE_TOB=${10}

PVALUE_BUILD_USAGE_LOG="${OUTPREFIX}.pvalue.build.usage.log"

### print all input arguments for debugging
echo "GT: $GT"
echo "ASTRAL_TREE: $ASTRAL_TREE"
echo "OUTPREFIX: $OUTPREFIX"
echo "TAXA: $TAXA"
echo "NET_ID: $NET_ID"
echo "TRUE_NET: $TRUE_NET"
echo "MODE: $MODE"
echo "ALPHA: $ALPHA"
echo "BETA: $BETA"
echo "TRUE_TOB: $TRUE_TOB"

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

########## compute annotation tree ##########

OUTPUT_PVALUE="${OUTPREFIX}.pvalue_build.tre"

if [[ "$MODE" == "3f1a" ]]; then
    TREE_QMC_ARGS="--3f1a"
else
    TREE_QMC_ARGS=""
fi

${TIMER} -v -o "${PVALUE_BUILD_USAGE_LOG}" \
    "${TREE_QMC}" \
    --blobsearchonly "${ASTRAL_TREE}" \
    --store_pvalue \
    ${TREE_QMC_ARGS} \
    --override \
    -i "${GT}" \
    -o "${OUTPUT_PVALUE}"

TIME_OF_AT_BUILDING=$(
    wall_time_seconds "${PVALUE_BUILD_USAGE_LOG}" | awk '{print $2}'
)

############ compute TOB ###############

OUTPUT_TOB="${OUTPREFIX}_alpha${ALPHA}_beta${BETA}.tob"
TOB_USAGE_LOG="${OUTPREFIX}_alpha${ALPHA}_beta${BETA}.tob.usage.log"

${TIMER} -v -o "${TOB_USAGE_LOG}" \
    "${TREE_QMC}" \
    --alpha "${ALPHA}" \
    --beta "${BETA}" \
    --load_pvalue \
    -i "${OUTPUT_PVALUE}" \
    -o "${OUTPUT_TOB}" \
    --override

TIME_OF_TOB_BUILDING=$(
    wall_time_seconds "${TOB_USAGE_LOG}" | awk '{print $2}'
)

TOTAL_TIME=$(awk \
    -v t1="$TIME_OF_AT_BUILDING" \
    -v t2="$TIME_OF_TOB_BUILDING" \
    'BEGIN{printf "%.3f", t1 + t2}')

############ compare TOB ###############

TREE_COMPATOR=/fs/cbcb-lab/ekmolloy/jdai123/network-cycle-sorting-study/tools/compare_two_tree.py

TOB_CSV="${OUTPUT_TOB}.csv"

COMP_LINE="$(
    python3 "${TREE_COMPATOR}" \
        -t1 "${TRUE_TOB}" \
        -t2 "${OUTPUT_TOB}"
)"

echo "taxa,network_id,true_tob,estimated_tob,method,alpha,beta,n1,n2,nl,i1,i2,fn,fp,fnr,fpr,(fnr+fpr)/2,wall_clock_seconds" > "${TOB_CSV}"

echo "${TAXA},${NET_ID},${TRUE_TOB},${OUTPUT_TOB},TOB-QMC-${MODE},${ALPHA},${BETA},${COMP_LINE},${TOTAL_TIME}" >> "${TOB_CSV}"
