#!/bin/bash
# run_pipeline.sh

# This is the main script for running the bulk RNAseq pipeline
# The setup.cfg and metadata.tsv files should be prepared prior to running this script

# script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#cd "$SCRIPT_DIR"

set -euo pipefail # Exit on error, undefined variable, or failed pipe

# Set up early tmp log file

# Backup original stdout/stderr
exec 3>&1 4>&2

# Create early log and redirect stdout/stderr to it, but also echo to original stdout
early_log=".preconfig_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee "$early_log" >&3) 2>&1

echo "Start Time: $(date)"
echo "Host: $(hostname)"
echo "User: $(whoami)"
echo "Job ID: ${PBS_JOBID:-N/A}"


#############
# Functions #
#############

confirm_continue() {
  $SCRIPT_DIR/scripts/confirm_continue.sh "$1" || exit 1
}

usage() {
  echo "Usage: $0 -c <config_file>" >&2
  exit 1
}

############################
# Handling input arguments #
############################

# Default values
config_file="setup.cfg"

# getopt setup
TEMP=$(getopt -o c: --long config: -n 'run_pipeline.sh' -- "$@") || usage

eval set -- "$TEMP"

while true; do
  case "$1" in
    -c|--config)
      config_file="$2"; shift 2;;
    --) shift; break;;
    *) usage;;
  esac
done

####################
# Read Config File #
####################

echo -e "\n\n---------------------------"
echo -e "Reading Config File"
echo -e "---------------------------\n\n"

if [ ! -f "$config_file" ] || [ ! -r "$config_file" ] || [ ! -s "$config_file" ]; then
  echo 'Config file does not exist, is not readable, or is empty' >&2
  exit 1
fi

source "$config_file"

# Default fallbacks
: "${threads:=4}"

if [[ -z "$output_dir" ]]; then
  echo "Error: output_dir not set in config file" >&2
  exit 1
fi

if [ "$interactive" = true ]; then

  if [ -d "$output_dir" ]; then
    echo "Output directory already exists and files may be overwritten."
    confirm_continue "Are you sure you want to continue?"
    mkdir -p "$output_dir"
  fi

fi

# Setup logging in output dir
mkdir -p "$output_dir/logs"
timestamp=$(date +%Y%m%d_%H%M%S)
final_log="$output_dir/logs/pipeline_${timestamp}.log"

# Restore original stdout/stderr before redirecting again
exec 1>&3 2>&4

# Start final logging (tee to both file and terminal)
exec > >(tee "$final_log") 2>&1

# Append early log contents to file only (not reprinted to screen)
cat "$early_log" >> "$final_log"
rm -f "$early_log"

#####################
# Setup Environment #
#####################

if [ "$user_environment" = "gadi" ]; then
  source env_setup/env_gadi.sh
fi

#######################
# Get FASTQ Filenames #
#######################

# First getting all filenames
filename_col=$(head -1 "$metadata" | tr '\t' '\n' | grep -nx 'filename' | cut -d: -f1)

# Read all filenames into an array
mapfile -t all_files < <(tail -n +2 "$metadata" | cut -f"$filename_col")


# Next need to get forward and reverse reads per sample
# Get header line from the metadata
header=$(head -n 1 "$metadata")

# Get column numbers by header name
sample_col=$(echo "$header" | tr '\t' '\n' | grep -nx 'sample_ID' | cut -d: -f1)
filename_col=$(echo "$header" | tr '\t' '\n' | grep -nx 'filename' | cut -d: -f1)

# Safety check: ensure required columns were found
if [[ -z "$sample_col" || -z "$filename_col" ]]; then
  echo "Error: sample_ID and/or filename column not found in metadata." >&2
  exit 1
fi


declare -A r1_files   # sample_ID → R1 filepath
declare -A r2_files   # sample_ID → R2 filepath

while IFS=$'\t' read -r line; do
  # Extract sample_ID and filename from this line using awk based on column number
  sample=$(echo "$line" | awk -v col=$sample_col -F'\t' '{print $col}')
  filename=$(echo "$line" | awk -v col=$filename_col -F'\t' '{print $col}')
  
  # Determine if this is an R1 or R2 read based on filename pattern
  if [[ "$filename" =~ _R1 ]]; then
    r1_files["$sample"]="$filename"
  elif [[ "$filename" =~ _R2 ]]; then
    r2_files["$sample"]="$filename"
  else
    echo "Warning: File for sample $sample is neither R1 nor R2: $filename" >&2
  fi
done < <(tail -n +2 "$metadata")

#-------------------------------------------

if [ "$run1_rawqc" = true ]; then

  ##############################
  # Step 1: Raw FASTQ QC Check #
  ##############################

  
  echo -e "\n\n---------------------------"
  echo "**Step 1: Raw read QC checking"
  echo -e "---------------------------\n\n"
  
  mkdir -p "$output_dir/qc/raw_reads"
  mkdir -p "$output_dir/logs/fastqc"
  
  # run fastqc on all samples in parallel
  $fastqc_path/fastqc \
  -t "$threads" \
  --outdir "$output_dir/qc/raw_reads" \
  "${all_files[@]}"

  $multiqc_path/multiqc "$output_dir/qc/raw_reads" -o "$output_dir/qc/raw_reads"
  
  
  echo "Raw read QC is available in $output_dir/qc/raw_reads"
  
  if [ "$interactive" = true ]; then
    echo
    echo "You have selected interactive mode, waiting for confirmation to continue..."
    confirm_continue "Do you want to continue?"
  fi
  
fi

#-------------------------------------------

if [ "$run2_trimming" = true ]; then
  #########################
  # Step 2: Trim Sequences #
  #########################
  
  
  echo -e "\n\n---------------------------"
  echo "**Step 2: Trimming adapter sequences and low quality reads"
  echo -e "---------------------------\n\n"
  
  echo "Using adapter sequences in $bbduk_adapters"
  echo "Using following parameters for trimming (can be changed in config file):"
  echo "  ktrim=$bbduk_ktrim"
  echo "  k=$bbduk_k"
  echo "  mink=$bbduk_mink"
  echo "  hdist=$bbduk_hdist"
  echo "  qtrim=$bbduk_qtrim"
  echo "  trimq=$bbduk_trimq"
  echo "  minlen=$bbduk_minlen"
  
  
  mkdir -p "$output_dir/intermediate_files"
  
  for sample in "${!r1_files[@]}"; do
  
    fwd="${r1_files[$sample]}"
    rev="${r2_files[$sample]}"
  
    echo "Trimming adapters in sample: $sample"
    echo "Forward read: $fwd"
    echo "Reverse read: $rev"
    
    $bbamp_path/bbduk.sh -Xmx20g \
      t=$threads \
      ref=$bbduk_adapters \
      in1=$fwd \
      in2=$rev \
      out1="$output_dir/intermediate_files/${sample}_R1_trim.fastq.gz" \
      out2="$output_dir/intermediate_files/${sample}_R2_trim.fastq.gz" \
      outm1="$output_dir/intermediate_files/${sample}_R1_trim_fail.fastq.gz" \
      outm2="$output_dir/intermediate_files/${sample}_R2_trim_fail.fastq.gz" \
      outs="$output_dir/intermediate_files/${sample}_pass_singletons.fastq.gz" \
      ktrim=$bbduk_ktrim \
      k=$bbduk_k \
      mink=$bbduk_mink \
      hdist=$bbduk_hdist \
      qtrim=$bbduk_qtrim \
      trimq=$bbduk_trimq \
      minlen=$bbduk_minlen \
      tpe tbo
      #tpe trims both paired reads to same length, tbo trims overlapping reads
  
  done
  
  echo "Trimmed fastq files are in $output_dir/intermediate_files"

fi

#-------------------------------------------

if [ "$run3_trimqc" = true ]; then

  ##################################
  # Step 3: Trimmed FASTQ QC Check #
  ##################################
  
  
  echo -e "\n\n---------------------------"
  echo "**Step 3: Trimmed read QC checking"
  echo -e "---------------------------\n\n"
  
  mkdir -p "$output_dir/qc/trimmed_reads"
  
  # collect all the trimmed FASTQ paths into an array
  trimmed_files=( "$output_dir"/intermediate_files/*_trim.fastq.gz )
  
  "$fastqc_path/fastqc" \
    -t "$threads" \
    --outdir "$output_dir/qc/trimmed_reads" \
    "${trimmed_files[@]}"
  
  #multiqc
  
  multiqc "$output_dir/qc/trimmed_reads" -o "$output_dir/qc/trimmed_reads"
  
  echo "Trimmed read QC is available in $output_dir/qc/raw_reads"
  
  if [ "$interactive" = true ]; then
    echo
    echo "You have selected interactive mode, waiting for confirmation to continue..."
    confirm_continue "Do you want to continue?"
  fi
  
fi

#-------------------------------------------



if [ "$run4_alignment" = true ]; then

  #####################
  # Step 4: Alignment #
  #####################
  
  
  echo -e "\n\n---------------------------"
  echo "**Step 4: Sequence Alignment"
  echo -e "---------------------------\n\n"
  
  
  
fi