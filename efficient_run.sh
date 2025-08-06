#efficient_run.sh
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Functions

usage() {
cat <<EOF >&2
Usage: $0 [OPTIONS]
  -c, --config      Full path to config file
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
    --pbs-walltime "3:00:00"
  
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
    --pbs-walltime "5:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_1:-}" ]]; then
      job_id_2=$(qsub -W depend=afterok:"$job_id_1" "$pbsdir/2_trimming.pbs")
    else
      job_id_2=$(qsub "$pbsdir/2_trimming.pbs")
    fi
  fi
fi
#---------------------------Step 3---------------------------#
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
    --pbs-walltime "3:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_2:-}" ]]; then
      job_id_3=$(qsub -W depend=afterok:"$job_id_2" "$pbsdir/3_trim_qc.pbs")
    else
      job_id_3=$(qsub "$pbsdir/3_trim_qc.pbs")
    fi
  fi
fi
#---------------------------Step 4o1---------------------------#
if [ "$run4o1_alignment_generate_ref" = true ]; then

  echo "Generating $pbsdir/4o1_generate_aligner_ref.pbs"
  
  star_ref_threads=$(calc_resource_per_job 12 1 "$threads")
  star_ref_mem=$(calc_resource_per_job 60 1 "$memory")
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps "4o1" \
    --pbs-script-name "$pbsdir/4o1_generate_aligner_ref.pbs" \
    --pbs-job-name "gen_ref" \
    --keep true \
    --threads "$star_ref_threads" \
    --pbs-mem "${star_ref_mem}GB" \
    --pbs-walltime "5:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_3:-}" ]]; then
      job_id_4o1=$(qsub -W depend=afterok:"$job_id_3" "$pbsdir/4o1_generate_aligner_ref.pbs")
    else
      job_id_4o1=$(qsub "$pbsdir/4o1_generate_aligner_ref.pbs")
    fi
  fi
fi
#---------------------------Step 4---------------------------#
if [ "$run4_alignment" = true ]; then
  echo "Generating $pbsdir/4_alignment.pbs"
  
  conc_jobs=$(calc_concurrent_jobs "$STAR_threads_per_job" $threads)
  star_threads=$(calc_resource_per_job "$STAR_threads_per_job" "$conc_jobs" "$threads")
  # smaller amount per process
  per_job_star_mem=$(calc_resource_per_job 5 "$conc_jobs" "$memory")
  # add larger armount for shared genome load based on index size
  index_dir="$STAR_index_dir/${read_lengths}bp"
  # measure number of GB for index
  base_index_gb=$(du -sBG --apparent-size "$index_dir" \
  | cut -f1 | sed 's/G$//')
               
  star_mem=$(( per_job_star_mem + base_index_gb + 5))
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps "4" \
    --pbs-script-name "$pbsdir/4_alignment.pbs" \
    --pbs-job-name "alignment" \
    --keep true \
    --threads "$star_threads" \
    --pbs-mem "${star_mem}GB" \
    --pbs-walltime "24:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous steps to finish
    if [[ -n "${job_id_3:-}" ]]; then
      job_id_4=$(qsub -W depend=afterok:"$job_id_3" "$pbsdir/4_alignment.pbs")
    elif [[ -n "${job_id_4o1}" ]]; then
      job_id_4=$(qsub -W depend=afterok:"$job_id_4o1" "$pbsdir/4_alignment.pbs")
    else
      job_id_4=$(qsub "$pbsdir/4_alignment.pbs")
    fi
  fi
fi
#---------------------------Step 5---------------------------#
if [ "$run5_alignqc" = true ]; then

  echo "Generating $pbsdir/5_align_qc.pbs"
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 5 \
    --pbs-script-name "$pbsdir/5_align_qc.pbs" \
    --pbs-job-name "alignQC" \
    --keep true \
    --threads 4 \
    --pbs-mem "8GB" \
    --pbs-walltime "2:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_4:-}" ]]; then
      job_id_5=$(qsub -W depend=afterok:"$job_id_4" "$pbsdir/5_align_qc.pbs")
    else
      job_id_5=$(qsub "$pbsdir/5_align_qc.pbs")
    fi
  fi
fi
#---------------------------Step 6o1---------------------------#
#if [ "$run6o1_rsem_generate_ref" = true ]; then

#fi
#---------------------------Step 6---------------------------#
if [ "$run6_quant" = true ]; then
  
  echo "Generating $pbsdir/6_quant.pbs"
  
  conc_jobs=$(calc_concurrent_jobs "$quant_threads_per_job" $threads)
  quant_threads=$(calc_resource_per_job "$quant_threads_per_job" "$conc_jobs" "$threads")
  
  if [ "$quantification" = "featureCounts" ]; then
    quant_mem=$(calc_resource_per_job 2 "$conc_jobs" "$memory")
  fi
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 6 \
    --pbs-script-name "$pbsdir/6_quant.pbs" \
    --pbs-job-name "quant" \
    --keep true \
    --threads $quant_threads \
    --pbs-mem "${quant_mem}GB" \
    --pbs-walltime "24:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_5:-}" ]]; then
      job_id_6=$(qsub -W depend=afterok:"$job_id_5" "$pbsdir/6_quant.pbs")
    else
      job_id_6=$(qsub "$pbsdir/6_quant.pbs")
    fi
  fi
  
fi
#---------------------------Step 7---------------------------#
if [ "$run7_consensusDE" = true ]; then
  echo "Generating $pbsdir/7_consensusDE.pbs"
  
  consensusDE_mem=$(calc_resource_per_job 1 "$n_samples" "$memory")
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 7 \
    --pbs-script-name "$pbsdir/7_consensusDE.pbs" \
    --pbs-job-name "consensusDE" \
    --keep true \
    --threads 12 \
    --pbs-mem "${consensusDE_mem}GB" \
    --pbs-walltime "24:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_6:-}" ]]; then
      job_id_7=$(qsub -W depend=afterok:"$job_id_6" "$pbsdir/7_consensusDE.pbs")
    else
      job_id_7=$(qsub "$pbsdir/7_consensusDE.pbs")
    fi
  fi
fi
#---------------------------Step 8---------------------------#
if [ "$run8_analyse_expression" = true ]; then
  echo "Generating $pbsdir/8_analyse_expression.pbs"
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 8 \
    --pbs-script-name "$pbsdir/8_analyse_expression.pbs" \
    --pbs-job-name "analyse_expression" \
    --keep true \
    --threads 4 \
    --pbs-mem "8GB" \
    --pbs-walltime "4:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_7:-}" ]]; then
      job_id_8=$(qsub -W depend=afterok:"$job_id_7" "$pbsdir/8_analyse_expression.pbs")
    else
      job_id_8=$(qsub "$pbsdir/8_analyse_expression.pbs")
    fi
  fi
fi
#---------------------------Step 9---------------------------#
if [ "$run9_enrichment_analysis" = true ]; then
  echo "Generating $pbsdir/9_enrichment_analysis.pbs"
  
  $SCRIPT_DIR/qsub_pipeline.sh \
    -c $config_file \
    --steps 9 \
    --pbs-script-name "$pbsdir/9_enrichment_analysis.pbs" \
    --pbs-job-name "enrichment" \
    --keep true \
    --threads 12 \
    --pbs-mem "16GB" \
    --pbs-walltime "8:00:00"
  
  if [ "$do_submit" == true ]; then
    # wait for previous step to finish
    if [[ -n "${job_id_8:-}" ]]; then
      job_id_9=$(qsub -W depend=afterok:"$job_id_8" "$pbsdir/9_enrichment_analysis.pbs")
    else
      job_id_9=$(qsub "$pbsdir/9_enrichment_analysis.pbs")
    fi
  fi
fi