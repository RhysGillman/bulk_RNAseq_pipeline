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


# setup resource monitoring log

RESOURCE_LOG="$output_dir/logs/resource_utilisation_${timestamp}.log"
: > "$RESOURCE_LOG"      # truncate or create



echo -e "\n\n"
echo "Config file read successfully"
echo "Running the following steps of pipeline"
echo "Raw Read QC: $run1_rawqc"
echo "Read Trimming: $run2_trimming"
echo "Trimmed Read QC: $run3_trimqc"
echo "Generate Aligner Index (optional): $run4o1_alignment_generate_ref"
echo "Alignment: $run4_alignment"
echo "Alignment QC: $run5_alignqc"
echo "Generate RSEM Index (Optional): $run6o1_rsem_generate_ref"
echo "Read Quantification: $run6_quant"
echo -e "\n\n"

##################
# More Functions #
##################

monitor_process() {
  local process_name="$1"; shift
  echo "===== $process_name START: $(date '+%F %T') =====" >> "$RESOURCE_LOG"
  
  # uses the time function to capture resource utilisation information
  # output is teed to both stdout and and the log file
  # grep filtered to only key information
  
  # create a tmp file for logging
  tmp_time=$(mktemp)

  # run the command under time function, send the command's stdout+stderr to the main log
  # capture only the time -v verbose stats into tmp $tmp_time via -o
  /usr/bin/time -v -o "$tmp_time" --append "$@"  > /dev/null 2>&1 || true

  # pull out relevant fields
  grep -E \
    'Command being timed:|User time \(seconds\):|System time \(seconds\):|Percent of CPU this job got:|Elapsed \(wall clock\) time|Maximum resident set size \(kbytes\):|File system inputs:|File system outputs:' \
    "$tmp_time" \
    >> "$RESOURCE_LOG" \
    || true
  
  # remove tmp log file
  rm "$tmp_time"
  
  echo "===== $process_name END:   $(date '+%F %T') =====" >> "$RESOURCE_LOG"
}



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

  
  echo -e "\n\n-------------------------------"
  echo "**Step 1: Raw read QC checking"
  echo -e "-------------------------------\n\n"
  
  mkdir -p "$output_dir/qc/raw_reads"
  
  # run fastqc on all samples in parallel
  monitor_process "Raw Read QC" \
    "$fastqc_path/fastqc" \
    -t "$threads" \
    --outdir "$output_dir/qc/raw_reads" \
    "${all_files[@]}"

  "$multiqc_path/multiqc" "$output_dir/qc/raw_reads" -o "$output_dir/qc/raw_reads"
  
  
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
  
  
  echo -e "\n\n--------------------------------------------------------------"
  echo "**Step 2: Trimming adapter sequences and low quality reads"
  echo -e "--------------------------------------------------------------\n\n"
  
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
    
    monitor_process "Sequence Trimming" \
      "$bbmap_path/bbduk.sh" -Xmx5g \
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
  
  
  echo -e "\n\n---------------------------------"
  echo "**Step 3: Trimmed read QC checking"
  echo -e "---------------------------------\n\n"
  
  mkdir -p "$output_dir/qc/trimmed_reads"
  
  # collect all the trimmed FASTQ paths into an array
  trimmed_files=( "$output_dir"/intermediate_files/*_trim.fastq.gz )
  
  monitor_process "Trimmed Read QC" \
    "$fastqc_path/fastqc" \
      -t "$threads" \
      --outdir "$output_dir/qc/trimmed_reads" \
      "${trimmed_files[@]}"
  
  #multiqc
  
  "$multiqc_path/multiqc" "$output_dir/qc/trimmed_reads" -o "$output_dir/qc/trimmed_reads"
  
  echo "Trimmed read QC is available in $output_dir/qc/raw_reads"
  
  if [ "$interactive" = true ]; then
    echo
    echo "You have selected interactive mode, waiting for confirmation to continue..."
    confirm_continue "Do you want to continue?"
  fi
  
fi

#-------------------------------------------

if [ "$run4o1_alignment_generate_ref" = true ]; then

  #############################################
  # Step 4(optional): Generate STAR Reference #
  #############################################
  
  if [ "$interactive" = true ]; then
  
    if [ -d "$STAR_index_dir" ]; then
      echo
      echo "You have chosen to create a STAR genome index, but the STAR index"
      echo "directory already exists"
      confirm_continue "Do you want to continue?"
    fi
    
  fi
  
  echo -e "\n\n----------------------------------------------------"
  echo "**Step 4 (Optional Step 1): STAR Index Generation"
  echo -e "----------------------------------------------------\n\n"
  
  
  mkdir -p "$STAR_index_dir/${read_lengths}bp"
  
  monitor_process "Generate STAR Index" \
    "$STAR_path/STAR" --runMode genomeGenerate \
      --genomeFastaFiles "$ref_fasta" \
      --sjdbGTFfile "$ref_gtf" \
      --genomeDir "$STAR_index_dir/${read_lengths}bp" \
      --runThreadN $threads \
      --sjdbOverhang $((read_lengths-1))
  
  echo "STAR index has been generated in $STAR_index_dir/${read_lengths}bp"
  
fi

#-------------------------------------------

if [ "$run4_alignment" = true ]; then

  #####################
  # Step 4: Alignment #
  #####################
  
  
  echo -e "\n\n---------------------------"
  echo "**Step 4: Sequence Alignment"
  echo -e "---------------------------\n\n"
  
  echo "Using STAR index in $STAR_index_dir/${read_lengths}bp"
  echo "Using following parameters for STAR alignment (can be changed in config file):"
  echo "sjdbOverhang=$((read_lengths-1))"
  echo "alignSJDBoverhangMin=$STAR_alignSJDBoverhangMin"
  echo "outFilterMultimapNmax=$STAR_outFilterMultimapNmax"
  echo "genomeSAindexNbases=$STAR_genomeSAindexNbases"
  echo -e "\n\n"
  
  
  # determine how many jobs to run in parallel, based on $threads
  MIN_THREADS_PER_JOB=$STAR_threads_per_job
  if [ "$threads" -lt "$MIN_THREADS_PER_JOB" ]; then
    STAR_jobs=1
    threads_per_job=$threads
  else
    STAR_jobs=$(( threads / MIN_THREADS_PER_JOB ))
    threads_per_job=$MIN_THREADS_PER_JOB
  fi
  echo "Running up to $STAR_jobs STAR jobs in parallel, each with $threads_per_job threads."

  # loop over samples and launch STAR in background, throttled to $star_jobs
  for sample in "${!r1_files[@]}"; do
    fwd="$output_dir/intermediate_files/${sample}_R1_trim.fastq.gz"
    rev="$output_dir/intermediate_files/${sample}_R2_trim.fastq.gz"
	  
	  echo "Beginning alignment for ${sample} using:"
	  echo "$fwd"
	  echo "$rev"
	  echo -e "\n\n"
    (
  	monitor_process "Read Alignment" \
      "$STAR_path/STAR" --genomeDir "$STAR_index_dir/${read_lengths}bp" \
        --readFilesIn "$fwd" "$rev" \
        --runThreadN $threads_per_job \
        --outSAMattrRGline ID:"$sample" SM:"$sample"_l1 PL:ILLUMINA \
        --sjdbOverhang $((read_lengths-1)) \
        --alignSJDBoverhangMin $STAR_alignSJDBoverhangMin \
        --outFilterMultimapNmax $STAR_outFilterMultimapNmax \
        --readFilesCommand zcat \
        --outFileNamePrefix "$output_dir/intermediate_files/${sample}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --outMultimapperOrder Random \
        --genomeSAindexNbases $STAR_genomeSAindexNbases
      
    echo -e "Completed STAR Alignment, basic info below, more in $output_dir/intermediate_files/${sample}_Log.final.out\n"    
    # print important log information
    grep -E "Uniquely mapped reads|reads mapped to multiple loci|reads unmapped" "$output_dir/intermediate_files/${sample}_Log.final.out"
    ) &
    # throttle: wait if too many RSEM jobs are running
    while [ "$(jobs -rp | wc -l)" -ge "$STAR_jobs" ]; do
      sleep 5
    done
    
    
    
  done

  wait  # ensure all background jobs finish before proceeding
  
  echo "Alignment files are available in $output_dir/intermediate_files/"

	for sample in "${!r1_files[@]}"; do
    samtools index "$output_dir/intermediate_files/${sample}_Aligned.sortedByCoord.out.bam"
  done
  
fi

#-------------------------------------------

if [ "$run5_alignqc" = true ]; then

  ######################
  # Step 5 Alignment QC#
  ######################
  
  
  echo -e "\n\n---------------------------"
  echo "**Step 5: Alignment QC checking"
  echo -e "---------------------------\n\n"
 
  mkdir -p "$output_dir/qc/alignment/"
  
  "$multiqc_path/multiqc" -o "$output_dir/qc/alignment/" $output_dir/intermediate_files/
  
  if [ "$interactive" = true ]; then
    echo
    echo "You have selected interactive mode, waiting for confirmation to continue..."
    confirm_continue "Do you want to continue?"
  fi
  
  
fi


#-------------------------------------------

if [ "$run6o1_rsem_generate_ref" = true ]; then

  #############################################
  # Step 6(optional): Generate RSEM Reference #
  #############################################
  
  #if [ -d "$RSEM_index_dir" ]; then
  #    echo
  #    echo "You have chosen to create RSEM reference files, but the directory"
  #    echo "$RSEM_index_dir already exists"
  #    confirm_continue "Do you want to continue?"
  #fi
  
  echo -e "\n\n----------------------------------------------------"
  echo "**Step 6 (Optional Step 1): RSEM Index Generation"
  echo -e "----------------------------------------------------\n\n"
  
  mkdir -p "$RSEM_index_dir"

  "$rsem_path/rsem-prepare-reference" \
    --quiet \
  	--gtf "$ref_gtf" \
  	"$ref_fasta" \
  	"$RSEM_index_dir"
  
  echo "Finished generating RSEM reference in $RSEM_index_dir"
  
fi

#-------------------------------------------

if [ "$run6_quant" = true ]; then
  ##########################
  # Step 6: Quantify Reads #
  ##########################
  
  echo -e "\n\n---------------------------"
  echo "**Step 6: Quantify Reads"
  echo -e "---------------------------\n\n"
  
  # auto-detect strandedness
  if [ "$stranded_reads" = "auto" ]; then
    echo "Auto-detecting strandedness..."
    source "$RSeQC_path/bin/activate"
    
    keys=( "${!r1_files[@]}" )
    first_sample="${keys[0]}"
    first_sample_bam="$output_dir/intermediate_files/${first_sample}_Aligned.sortedByCoord.out.bam"
    s=$(./scripts/guess_strandedness.sh "$first_sample_bam" "$ref_exon_bed")
    if [ -z "$s" ]; then
      echo "Warning: strandedness auto‐detection returned empty, defaulting to unstranded" >&2
      s="unstranded"
    else
      echo "Detected strandedness as: $s"
    fi
    deactivate
  else
    s=$stranded_reads
  fi
  
  case "$s" in
    unstranded)         rsem_s="none"    ;;
    forward-stranded)   rsem_s="forward" ;;
    reverse-stranded)   rsem_s="reverse" ;;
    *) 
      echo "Warning: unknown strandedness '$s', defaulting to unstranded" >&2
      rsem_s="none"
      ;;
  esac

  # determine how many RSEM jobs to run in parallel, based on $threads
  MIN_THREADS_PER_JOB=$RSEM_threads_per_job
  if [ "$threads" -lt "$MIN_THREADS_PER_JOB" ]; then
    rsem_jobs=1
    threads_per_job=$threads
  else
    rsem_jobs=$(( threads / MIN_THREADS_PER_JOB ))
    threads_per_job=$MIN_THREADS_PER_JOB
  fi
  echo "Running up to $rsem_jobs RSEM jobs in parallel, each with $threads_per_job threads."

  # loop over samples and launch RSEM in background, throttled to $rsem_jobs
  for sample in "${!r1_files[@]}"; do
    echo -e "\nQuantifying reads for sample: $sample using RSEM\n"
    (
      monitor_process "RSEM Quantification: $sample" \
        "$rsem_path/rsem-calculate-expression" \
          -p "$threads_per_job" \
          --paired-end \
          --quiet \
          --alignments \
          --strandedness "$rsem_s" \
          --no-bam-output \
          "$output_dir/intermediate_files/${sample}_Aligned.toTranscriptome.out.bam" \
          "$RSEM_index_dir" \
          "$output_dir/intermediate_files/${sample}"
      
      echo -e "Completed RSEM for sample: $sample"
      
    ) &
    # throttle: wait if too many RSEM jobs are running
    while [ "$(jobs -rp | wc -l)" -ge "$rsem_jobs" ]; do
      sleep 5
    done
  done

  wait  # ensure all background RSEM jobs finish before proceeding
  
  echo "RSEM results files are available in $output_dir/intermediate_files/"
  
  
  case "$s" in
    unstranded)         featurecounts_s=0 ;;
    forward-stranded)   featurecounts_s=1 ;;
    reverse-stranded)   featurecounts_s=2 ;;
    *) 
      echo "Warning: unknown strandedness '$s', defaulting to unstranded" >&2
      rsem_s="none"
      ;;
  esac
  
  for sample in "${!r1_files[@]}"; do
    echo -e "\nQuantifying reads for sample: $sample using featureCounts\n"
    (
      monitor_process "FeatureCounts Quantification: $sample" \
          "$subread_path/featureCounts" \
          -p \
          -T "$threads_per_job" \
          -s $featurecounts_s \
          -t exon \
          -g gene_id \
          -a "$ref_gtf" \
          -o "$output_dir/intermediate_files/${sample}_featureCounts.output" \
          "$output_dir/intermediate_files/${sample}_Aligned.sortedByCoord.out.bam"
      echo -e "Completed featureCounts for sample: $sample"
    ) &
    
    
    # throttle: wait if too many featurecounts jobs are running
    while [ "$(jobs -rp | wc -l)" -ge "$rsem_jobs" ]; do
      sleep 5
    done
  done

  wait  # ensure all background featurecounts jobs finish before proceeding
  
  echo "featureCounts results files are available in $output_dir/intermediate_files/"
  
  
  for sample in "${!r1_files[@]}"; do
    awk 'NR>2 {print $1, $NF}' "$output_dir/intermediate_files/${sample}_featureCounts.output" > "$output_dir/intermediate_files/${sample}_featureCounts.txt"
  done
  
  
  
fi

#-------------------------------------------

if [ "$run7_consensusDE" = true ]; then

  #######################
  # Step 7: ConsensusDE #
  #######################
  
  echo -e "\n\n----------------------------------------"
  echo "**Step 7: Find DEGs using ConsensusDE"
  echo -e "----------------------------------------\n\n"
  
  echo "Summarising isoform-level counts to gene-level counts"
  
  mkdir -p "$output_dir/consensusDE/input"
  
  Rscript --no-save --no-restore "scripts/summarise_gene_counts.R" -i "$output_dir/intermediate_files" -o "$output_dir/consensusDE/input"
  
  
  
  
fi




