#efficient_run.sh
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Functions

usage() {
cat <<EOF >&2
Usage: $0 [OPTIONS]
  -c, --config      Path to config file (default: setup.cfg)
  -t, --threads     Total threads available for job
  -m, --memory      Total memory available for job (in GB)
  -p, --pbsdir      Directory to store PBS submission scripts in
  --submit          Directly submit all jobs
  --steps           1,2,3...     Comma-separated list of steps to run
EOF
exit 1
}


calc_resource_per_job() {
  local per_job=$1
  local n_jobs=$2
  local avail=$3
  local need=$(( per_job * n_jobs ))
  if (( need <= avail )); then
    echo "$need"
  else
    echo "$avail"
  fi
}

calc_concurrent_jobs() {
  local per_job=$1
  local avail=$2
  local conc_jobs=$(( avail / per_job ))
  if (( conc_jobs > 0 )); then
    echo "$conc_jobs"
  else
    echo 1
  fi
}

############################
# Handling input arguments #
############################

# Default values
config_file="setup.cfg"
pbsdir="./"
do_submit=false
cli_steps="1,2,3,4,4o1,5,6o1,6,7,8,9"

# getopt setup
TEMP=$(getopt -o c:,t:,m:,p: \
  --long config:,threads:,memory:,pbsdir:,steps:,submit \
  -n 'efficient_run.sh' -- "$@") || usage
eval set -- "$TEMP"

while true; do
  case "$1" in
    -c|--config)      config_file="$2";       shift 2 ;;
    -t|--threads)     threads="$2";           shift 2 ;;
    -m|--memory)      memory="$2";            shift 2 ;;
    -p|--pbsdir)      pbsdir="$2";            shift 2 ;;
    --submit)         do_submit=true;         shift 1 ;;
    --steps)          cli_steps="$2";         shift 2 ;;
    --)               shift; break ;;
    *)                usage ;;
  esac
done

source "$config_file"

#echo "DEBUG: threads=$threads  memory=$memory  pbsdir=$pbsdir  steps=$cli_steps  submit=$do_submit"

mkdir -p $pbsdir

# Handling CLI selected steps to run

if [ -n "$cli_steps" ]; then
  # First, turn *all* main steps off
  run1_rawqc=false
  run2_trimming=false
  run3_trimqc=false
  run4_alignment=false
  run4o1_alignment_generate_ref=false
  run5_alignqc=false
  run6_quant=false
  run6o1_rsem_generate_ref=false
  run7_consensusDE=false
  run8_analyse_expression=false
  run9_enrichment_analysis=false

  # Parse comma-separated list
  IFS=',' read -ra STEP_ARR <<< "$cli_steps"
  for st in "${STEP_ARR[@]}"; do
    case "$st" in
      1) run1_rawqc=true                              ;;
      2) run2_trimming=true                           ;;
      3) run3_trimqc=true                             ;;
      4) run4_alignment=true                          ;;
      4o1) run4o1_alignment_generate_ref=true         ;;
      5) run5_alignqc=true                            ;;
      6) run6_quant=true                              ;;
      6o1) run6o1_rsem_generate_ref=true              ;;
      7) run7_consensusDE=true                        ;;
      8) run8_analyse_expression=true                 ;;
      9) run9_enrichment_analysis=true                ;;
      *) echo "Error: unknown step '$st'" >&2; exit 1 ;;
    esac
  done
fi


# Work Out Sample Size

header=$(head -n 1 "$metadata")
# Get column numbers by header name
sample_col=$(echo "$header" | tr '\t' '\n' | grep -nx 'sample_ID' | cut -d: -f1)
filename_col=$(echo "$header" | tr '\t' '\n' | grep -nx 'filename' | cut -d: -f1)

# Safety check: ensure required columns were found
if [[ -z "$sample_col" || -z "$filename_col" ]]; then
  echo "Error: sample_ID and/or filename column not found in metadata." >&2
  exit 1
fi

n_files=$(tail -n +2 "$metadata" | cut -f"$filename_col" | sort | uniq | wc -l)
n_samples=$(tail -n +2 "$metadata" | cut -f"$sample_col" | sort | uniq | wc -l)

echo "Detected $n_files files from metadata"
echo "Detected $n_samples samples from metadata"

#---------------------------Step 1---------------------------#
if [ "$run1_rawqc" = true ]; then
  conc_jobs=$(calc_concurrent_jobs 1 $threads)
  fastqc_threads=$(calc_resource_per_job 1 "$conc_jobs" "$threads")
  fastqc_mem=$(calc_resource_per_job 1 "$conc_jobs" "$memory")
  
  echo "Generating $pbsdir/1_raw_fastqc.pbs"
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 1 \
    --pbs-script-name "$pbsdir/1_raw_fastqc.pbs" \
    --pbs-job-name "fastqc" \
    --keep true \
    --threads "$fastqc_threads" \
    --pbs-mem "${fastqc_mem}GB" \
    --pbs-walltime "1:00:00"
  
  if [ "$do_submit" == true ]; then
    job_id_1=$(qsub "$pbsdir/1_raw_fastqc.pbs")
  fi
fi
#---------------------------Step 2---------------------------#
if [ "$run2_trimming" = true ]; then
  conc_jobs=$(calc_concurrent_jobs "$bbduk_threads_per_job" $threads)
  bbduk_threads=$(calc_resource_per_job "$bbduk_threads_per_job" "$conc_jobs" "$threads")
  bbduk_mem=$(calc_resource_per_job 2 "$conc_jobs" "$memory")
  
  echo "Generating $pbsdir/2_trimming.pbs"
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 2 \
    --pbs-script-name "$pbsdir/2_trimming.pbs" \
    --pbs-job-name "trim" \
    --keep true \
    --threads "$bbduk_threads" \
    --pbs-mem "${bbduk_mem}GB" \
    --pbs-walltime "4:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_1:-}" ]]; then
      job_id_2=$(qsub -W depend=afterok:"$job_id_1" "$pbsdir/2_trimming.pbs")
    else
      job_id_2=$(qsub "$pbsdir/2_trimming.pbs")
    fi
  fi
fi

if [ "$run3_trimqc" = true ]; then
  conc_jobs=$(calc_concurrent_jobs 1 $threads)
  fastqc_threads=$(calc_resource_per_job 1 "$conc_jobs" "$threads")
  fastqc_mem=$(calc_resource_per_job 1 "$conc_jobs" "$memory")
  
  echo "Generating $pbsdir/3_trim_qc.pbs"
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 3 \
    --pbs-script-name "$pbsdir/3_trim_qc.pbs" \
    --pbs-job-name "trimQC" \
    --keep true \
    --threads "$fastqc_threads" \
    --pbs-mem "${fastqc_mem}GB" \
    --pbs-walltime "1:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_2:-}" ]]; then
      job_id_3=$(qsub -W depend=afterok:"job_id_2" "$pbsdir/3_trim_qc.pbs")
    else
      job_id_3=$(qsub "$pbsdir/3_trim_qc.pbs")
    fi
  fi
fi

#if [ "$run4o1_alignment_generate_ref" = true ]; then

#fi

#if [ "$run5_alignqc" = true ]; then

#fi

#if [ "$run6o1_rsem_generate_ref" = true ]; then

#fi

#if [ "$run6_quant" = true ]; then

#fi

#if [ "$run7_consensusDE" = true ]; then

#fi

#if [ "$run8_analyse_expression" = true ]; then

#fi

#if [ "$run9_enrichment_analysis" = true ]; then

#fi


# Parallel Steps
## BBduk, STAR, quant, 


# Simple Steps
# raw QC, trim QC, align QC, 