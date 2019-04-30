#!/bin/bash

###################################################################################################
# ENSURE THE PIPELINE IS CALLED CORRECTLY
###################################################################################################

# Set up a usage statement in case this program is called incorrectly
usage() { echo -e "\nERROR: Missing SRA accessions and/or input query and/or query type. \n\n" \
              "Make sure to provide one (or more) SRA run numbers separated by commas, \n" \
              "as well as a virus query (in fasta format), and indicate the query type \n" \
              "as either 'nucl' or 'prot' (do not include the quotes). \n\n" \
              "Proper usage: \n" \
              "$0 -s SRR10001,SRR10002,SRR... -q VIRUS_QUERY -t nucl|prot \n\n" \
              "Optional parameters: \n" \
                "-e (evalue, e.g. 100, 1, or 1e-99; [default = 1e-9]) \n" \
                "-m (maximum amount of memory to use [in GB]; [default=16] ) \n" \
                "-d (sets nucloetide program to discontiguous-megablast; [default=megablast] ) \n" \
                "-n (sets nucleotide program to blastn; [default=megablast] ) \n\n" \
              "Example of a complex run: \n" \
              "$0 -s SRX193147,SRX193148,SRX193149 -q tvv2_nt.fasta -t nucl -e 1e-3 -m 30 -d \n\n" \
              "Exiting program. Please retry with corrected parameters..." >&2; exit 1; }

# Make sure the pipeline is invoked correctly, with project and sample names
while getopts "s:q:t:e:m:dn" arg; do
        case ${arg} in
                s ) # Take in the sample name(s)
                  set -f
                  IFS=","
                  ALL_SAMPLES=(${OPTARG}) # call this when you want every individual sample
                        ;;
                q ) # Take in the name of the virus query file (fasta format)
                  VIRUS_QUERY=${OPTARG}
                        ;;
                t ) # Take in the query type (nucleotide or protein)
                  QUERY_TYPE=${OPTARG}
                        ;;
                e) # set evalue
                  E_VALUE=${OPTARG}
                        ;;
                m) # set max memory to use (in GB; if any letters are entered, discard those)
                  MEMORY_ENTERED=${OPTARG}
                  MEMORY_TO_USE=$(echo $MEMORY_ENTERED | sed 's/[^0-9]*//g')
                        ;;
                d) # switch to discontiguous_megablast
                  BLAST_TASK="dc-megablast"
                        ;;
                n) # switch to blastn
                  BLAST_TASK="blastn"
                        ;;
                * ) # Display help
                  usage
                        ;;
        esac
done
shift $(( OPTIND-1 ))
###################################################################################################

###################################################################################################
# PROCESS THE USER PROVIDED PARAMETERS
###################################################################################################
# If the pipeline is not called correctly, tell that to the user and exit
if [[ -z "${VIRUS_QUERY}" ]] || [[ -z "${ALL_SAMPLES}" ]] || [[ -z "${QUERY_TYPE}" ]] ; then
        usage
fi

# Retrieve the name of last sample (using older but cross-platform compatible BASH notation)
LAST_SAMPLE=${ALL_SAMPLES[${#ALL_SAMPLES[@]}-1]}

# Create a variable that other parts of this pipeline can use mostly for naming
SAMPLES="${ALL_SAMPLES[0]}-${LAST_SAMPLE}"

# Reset global expansion
set +f

# Handle the query type provided by the user, using that to determine which type of blast to use
## Nucleotide query
if [[ ${QUERY_TYPE} == 'nucl' ]]; then
    BLAST_TYPE='blastn'

    # If -d (dc_megablst) or -n (blastn) flags were no given by user, default to megablast
    if [[ -z "${BLAST_TASK}" ]]; then
        BLAST_TASK='megablast'
    fi

## Protein query
elif [[ ${QUERY_TYPE} == 'prot' ]]; then
    BLAST_TYPE='tblastn'
    BLAST_TASK='tblastn'

## Other/Error
else
        echo "QUERY_TYPE is ${QUERY_TYPE}"
        echo "QUERY_TYPE must be 'nucl' or 'prot' (do not include quotes)"
        usage
fi

# If e-value wasn provided by user, then set it to 1e-9
if [[ -z ${E_VALUE} ]]; then
    E_VALUE="1e-9"
fi
###################################################################################################

###################################################################################################
# ENSURE THAT ALL REQUIRED SOFTWARE IS INSTALLED
###################################################################################################
command -v seqtk > /dev/null || \
{ echo -e "This program requires 'seqtk'. \n" \
        "Please install from github.com/lh3/seqtk and retry. \n" \
        "Exiting..."; exit 2
}

{ command -v makeblastdb > /dev/null && \
  command -v ${BLAST_TYPE} > /dev/null
} \
|| \
{ echo -e "This program requires 'makeblastdb'. \n" \
        "Please install from NCBI website or NCBI github and retry. \n" \
        "Exiting..."; exit 3
}
###################################################################################################

###################################################################################################
# Set up number of CPUs to use (use all available) and RAM (pls configure)
###################################################################################################
# Use `nproc` if installed (Linux or MacOS with gnu-core-utils); otherwise use `systctl`
{ \
    command -v nproc > /dev/null && \
    NUM_THREADS=`nproc` && \
    echo "Number of processors available (according to nproc): ${NUM_THREADS}"; \
} || \
{ \
    command -v sysctl > /dev/null && \
    NUM_THREADS=`sysctl -n hw.ncpu` && \
    echo "Number of processors available (according to sysctl): ${NUM_THREADS}";
}

# Set memory usage
if [[ -z ${MEMORY_TO_USE} ]]; then
    echo "No memory limit set by user. Defaulting to 16GB"
    MEMORY_TO_USE="16"
fi
###################################################################################################

###################################################################################################
# CREATE DIRECTORIES AND PREPARE NAMES FOR BLAST
###################################################################################################

# Create a directory to run & store the BLAST files
mkdir -p ${SAMPLES}

# Create names for BLAST output file:
## truncates file path, leaving just the filename itself
VIRUS_QUERY_FILE=${VIRUS_QUERY##*/}

## eliminates file extension, giving a cleaner name for blast
BLAST_NAME_VIRUS_QUERY=${VIRUS_QUERY_FILE%.*}

# Create log file
readonly LOG_FILE="${SAMPLES}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.log"
touch ${LOG_FILE}

# Copy initial launch command into the log
echo -e "sab was launched with the following command: \n $0 $@ \n" > ${LOG_FILE}

# Set directory to save SRA files
SRA_DIR="${HOME}/Documents/Research/sra/"
mkdir -p ${SRA_DIR}

# Make directory to save resulting BLASTDB
BLAST_DB_DIR=${SRA_DIR}/blastdbs
mkdir -p ${BLAST_DB_DIR}

# Read inputs back to the user and store them in the log
echo -e "\n" \
        "SRA Accessions provided: ${ALL_SAMPLES[@]} \n" \
        "Virus query file provided: ${VIRUS_QUERY} \n" \
        "Molecule type (nucl or prot) of input query: ${QUERY_TYPE} \n" \
        "e-value: ${E_VALUE} \n" \
        "Blast program: ${BLAST_TYPE} > ${BLAST_TASK} \n" \
        "Number of processors to use: ${NUM_THREADS} \n" \
        "Memory limit: ${MEMORY_TO_USE}GB \n\n"| tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# Download SRA files
###################################################################################################

# Protect against from empty/unset variables (waited to now bc of all the parameter setting)
set -u

# Add the download from SRA step to the timelog file
echo -e "Downloading input FASTQs from the SRA at: `date` \n" | tee -a ${LOG_FILE}

# Download fastq files from the SRA
for SAMPLE in ${ALL_SAMPLES[@]};
  do \
     echo -e "Downloading FASTQ file(s) for ${SAMPLE} at: \n `date`" | tee -a ${LOG_FILE};
     fasterq-dump \
     --split-3 \
     -t /tmp \
     --progress --verbose \
     --skip-technical --rowid-as-name --print-read-nr \
     --threads=${NUM_THREADS} \
     --mem=${MEMORY_TO_USE}"GB" \
     --bufsize=1000MB \
     --curcache=1000MB \
     --outdir ${SRA_DIR} \
     ${SAMPLE} >> ${LOG_FILE};
  done

echo -e "\n Finished downloading FASTQs from the SRA at: \n" \
        "`date` \n" | tee -a ${LOG_FILE}

echo -e "These SRA FASTQ files are located at: ${SRA_DIR} \n\n" | tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# MAKE BLASTDB FROM SRA FILES
###################################################################################################
# Devise name for blast database from every SRA given (eg: SRA1_SRA2_SRA3_db)
BLAST_DB_NAME="`echo -e ${ALL_SAMPLES[@]} | tr ' ' '_'`_db"

# Concatenate all SRA FASTQs together to make blastdb from & extract hits from
CONCATENATED_FASTQ=${BLAST_DB_DIR}/${BLAST_DB_NAME}.fq
CONCATENATED_FASTA=${BLAST_DB_DIR}/${BLAST_DB_NAME}.fasta

# Check to see if an intact blastdb & concatenated FASTA exist; if so, skip blastdb building step
blastdbcmd -db ${BLAST_DB_DIR}/${BLAST_DB_NAME} -info >> ${LOG_FILE} 2>&1
if [[ $? -eq 0 && -f ${CONCATENATED_FASTA} ]] ; then
    echo -e "BLAST DB for these SRA accessions already exists. \n" \
            "Continuing to BLAST alignment..." | tee ${LOG_FILE}

else
  # Put starting time of blastdb building in log file
  echo -e "Building BLAST database from SRA files at: \n `date`" | tee -a ${LOG_FILE}

  # If there is already a concat fastq, then it could be a partial file, so delete it
  if [[ -f ${CONCATENATED_FASTQ} ]]; then rm ${CONCATENATED_FASTQ}; fi

  # For each SRA accession, add its reads to the master concatendated fastq file
  for SRA_FASTQ in ${ALL_SAMPLES[@]}
    do \
      cat ${SRA_DIR}/${SRA_FASTQ}*fastq >> ${CONCATENATED_FASTQ};
    done

  # Transform that FASTQ into FASTA so it's readable by BLAST
  seqtk seq -A ${CONCATENATED_FASTQ} > ${CONCATENATED_FASTA}

  # Create the blast databse
  makeblastdb \
  -dbtype nucl \
  -in  ${CONCATENATED_FASTA} \
  -title ${BLAST_DB_NAME} \
  -out ${BLAST_DB_DIR}/${BLAST_DB_NAME} \
  -logfile ${BLAST_DB_DIR}/${SAMPLES}_makeblastdb.log

  # Delete temporary concatenated fastq object
  rm ${CONCATENATED_FASTQ}

  # Put finishing time of blastdb building in log file
  echo -e "Completed buiding BLAST database at: `date` \n" \
          "Log file for this step available at: ${BLAST_DB_DIR}/${SAMPLES}_makeblastdb.log" |
          tee -a ${LOG_FILE}
fi
###################################################################################################

###################################################################################################
# RUN BLAST
###################################################################################################
# Print time that blast starts and write to log file
echo -e "Began running ${BLAST_TYPE} with samples ${SAMPLES} at:" | tee -a ${LOG_FILE}
date | tee -a ${LOG_FILE}

# Run blast
${BLAST_TYPE} \
-task ${BLAST_TASK} \
-db ${BLAST_DB_DIR}/${BLAST_DB_NAME} \
-query ${VIRUS_QUERY} \
-out ${SAMPLES}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.results.txt \
-outfmt "6 qseqid sseqid evalue" \
-num_threads ${NUM_THREADS} \
-evalue ${E_VALUE} \
-max_target_seqs 100000000

# Print time completed and write to log file as well
echo -e "Finished running ${BLAST_TYPE} at:" | tee -a ${LOG_FILE}
date | tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# Create FASTA sequence file of hits
###################################################################################################
seqtk subseq  ${CONCATENATED_FASTA} \
  <(cut -f 2 ${SAMPLES}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.results.txt | sort -u) > \
  ${SAMPLES}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.fasta

# Print number of hits and save to log file
echo -e "\nNumber of hits, saved in fasta file:" | tee -a ${LOG_FILE}
grep -c "^>" ${SAMPLES}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.fasta | tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# STAND-ALONE-BLAST (sab) FINISHED SUCCESSFULLY
###################################################################################################
echo -e "\n StandAloneBlast (sab) finished successfully at: \n" \
        "`date`" | tee -a ${LOG_FILE}
