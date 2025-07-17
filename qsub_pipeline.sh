#!/usr/bin/env bash

### Defaults for PBS options ###
pbs_project="pq84"    # --pbs-project
pbs_queue="normalbw"  # --pbs-queue
walltime="24:00:00"   # --walltime
mem="64GB"            # --mem
jobfs="5GB"           # --jobfs
storage="gdata/u86+scratch/u86+gdata/pq84+scratch/pq84"  # --storage
mail_opts="abe"       # --mail-opts
job_name="RNAseq_pipeline" # --job-name
mail_user="rhys.gillman@jcu.edu.au"          # --mail-user

###  Pull out ONLY PBS-flags ###
run_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    # PBS flags (consume but don’t forward)
    --pbs-project)   pbs_project="$2";   shift 2 ;;
    --pbs-queue)     pbs_queue="$2";     shift 2 ;;
    --pbs-walltime)  walltime="$2";      shift 2 ;;
    --pbs-mem)       mem="$2";           shift 2 ;;
    --pbs-jobfs)     jobfs="$2";         shift 2 ;;
    --pbs-storage)   storage="$2";       shift 2 ;;
    --pbs-mail-opts) mail_opts="$2";     shift 2 ;;
    --pbs-job-name)  job_name="$2";      shift 2 ;;
    --pbs-mail-user) mail_user="$2";     shift 2 ;;
    # anything else is for the pipeline
    *) run_args+=("$1"); shift ;;
  esac
done

# force non-interactive
run_args+=(--interactive false)

### Extract ncpus from run_args (the --threads value) ###
ncpus=""
for ((i=0;i<${#run_args[@]};i++)); do
  if [[ "${run_args[i]}" == "--threads" ]]; then
    ncpus="${run_args[i+1]}"
    break
  fi
done
if [[ -z "$ncpus" ]]; then
  echo "Error: --threads is required" >&2
  exit 1
fi

### Generating PBS Submission Script ###

pbs_file=$(mktemp qsub_pipeline.XXXXXXXX.pbs)
cat > "$pbs_file" <<EOF
#!/usr/bin/env bash
#PBS -P $pbs_project
#PBS -q $pbs_queue
#PBS -l walltime=$walltime
#PBS -l mem=$mem
#PBS -l jobfs=$jobfs
#PBS -l ncpus=$ncpus
#PBS -l storage=$storage
#PBS -j oe
#PBS -m $mail_opts
#PBS -N $job_name
EOF

if [[ -n "$mail_user" ]]; then
  echo "#PBS -M $mail_user" >> "$pbs_file"
fi

cat >> "$pbs_file" <<EOF
shopt -s expand_aliases
source /etc/profile.d/modules.sh

echo "Job identifier is \$PBS_JOBID"
echo "Working directory is \$PBS_O_WORKDIR"

cd "\$PBS_O_WORKDIR"
bash "\$(pwd)/run_pipeline.sh" ${run_args[@]}
EOF

echo "Submitting pipeline job via qsub…"
jobid=$(qsub "$pbs_file")
rm -f "$pbs_file"
echo "Submitted job $jobid"