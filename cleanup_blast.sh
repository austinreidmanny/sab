#!/bin/bash

###################################################################################################
# ENSURE THE SCRIPT IS CALLED CORRECTLY
###################################################################################################
set -eo pipefail

usage() { echo -e "The objective of this script is to clean up the results of stand-alone-blast (sab). \n\n" \
                  "This program takes in a FASTA (given as output by sab) of reads that mapped \n" \
                  "to the previous reference and searches those again, this time against a \n" \
                  "reference nucleotide blast database. The output of this script is a \n" \
                  "tab-delimited text file with the following fields: \n" \
                  "(1) name of the read (2) title of best hit in NCBI (3) e-value score of match \n\n" \
                  "Usage: $0 -i input.fasta [options] \n\n" \
                  "Output: input.cleanup.results.txt \n\n" \
                  "Optional parameters: \n" \
                        "-e (evalue, e.g. 100, 1, or 1e-99; [default = 10]) \n" \
                        "-m (maximum amount of memory to use [in GB]; [default=16] ) \n" \
                        "-p (path to directory for database; " \
                            "[default='~/Documents/Research/sra/blastdbs/tvv_db'] ) \n" \
                        "-n (sets nucleotide program to blastn; [default= dc-megablast] ) \n" \
                        "-g (sets nucleotide program to megablast; [default= dc-megablast] ) \n\n" \
                      "Example of a complex run: \n" \
                      "$0 -i input.fasta -e 1e-3 -m 26 -g \n\n" \
                      "Exiting program. Please retry with corrected parameters..." >&2 && exit 1
        }
###################################################################################################

###################################################################################################
# Note about the RefSeq Viral nucleotide database used
###################################################################################################
# User will need to download this database from ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ ,
# download the two files 'viral.1.1.genomic.fna.gz' and 'viral.2.1.genomic.genomic.fna.gz',
# unzip and concatenate them together into a single viral refseq fasta file, and use that as
# input into the NCBI `makeblastdb` tool to make a blast nucleotide database.
# For posterity, this tool was tested with these steps performed at 10:30AM on 2019-05-03
###################################################################################################

###################################################################################################
# TAKE IN THE USER-PROVIDED PARAMETERS
###################################################################################################
# Store all parameters so I can save them to the log later
ALL_PARAMETERS=$@

# Read in the user-provided parameters
while getopts "i:e:m:p:ng*:" arg; do
        case ${arg} in
                i ) # Take in the input fasta
                  SEQDUMP=${OPTARG}
                        ;;
                e ) # set evalue
                  E_VALUE=${OPTARG}
                        ;;
                m ) # set max memory to use (in GB; if any letters are entered, discard those)
                  MEMORY_ENTERED=${OPTARG}
                  MEMORY_TO_USE=$(echo $MEMORY_ENTERED | sed 's/[^0-9]*//g')
                        ;;
                p ) #set path to NCBI nt database
                  PATH_TO_NT_DB=${OPTARG}
                        ;;
                n ) # switch to blastn
                  BLAST_TASK="blastn"
                        ;;
                g ) # switch to megablast
                  BLAST_TASK="megablast"
                        ;;
                * ) # Display help
                  usage
                        ;;
        esac
done
shift $(( OPTIND-1 ))
###################################################################################################

###################################################################################################
# PROCESS THOSE PARAMETERS
###################################################################################################
# If no input is provided, tell that to the user and exit
if [[ -z "${SEQDUMP}" ]] ; then
    echo -e "\nERROR: No input detected. \n"
    usage
fi

# If e-value wasn't provided by user, then set it to 1e-9
if [[ -z ${E_VALUE} ]]; then
    E_VALUE="10"
fi

# If -n (blastn) or -g (megablast) flags were not given by user, default to dc-megablast
if [[ -z "${BLAST_TASK}" ]]; then
    BLAST_TASK="dc-megablast"
fi

# If path to NCBI nt database was not given, give default path
if [[ -z "${PATH_TO_NT_DB}" ]]; then
    PATH_TO_NT_DB="${HOME}/Documents/Research/sra/blastdbs/tvv_db"
fi
###################################################################################################

###################################################################################################
# Set up number of CPUs to use (use all available) and RAM (can be set by user, defaults to 16GB)
###################################################################################################
# CPUs (aka threads aka processors aka cores):
#   If 8 CPUs were used, BLAST fails & gives Segmentation Fault. Error stopped if <= 4 CPUs are used
#   Strategy: Use up to 4 CPUs, or maximum available if less than 4 CPUs available

# Use `nproc` if installed (Linux or MacOS with gnu-core-utils); otherwise use `systctl`
{ \
    command -v nproc > /dev/null && \
    MAX_NUM_THREADS=`nproc`
} || \
{ \
    command -v sysctl > /dev/null && \
    MAX_NUM_THREADS=`sysctl -n hw.ncpu`
}

# If maximum available threads is less than or equal to 4, use all threads; else use 4
if (( ${MAX_NUM_THREADS} > 4 )); then
    NUM_THREADS=4
elif (( ${MAX_NUM_THREADS} <= 4 )); then
    NUM_THREADS=${MAX_NUM_THREADS}
else
    echo "Error. Could not determine number of CPUs to use. Exiting..."
    exit 4
fi

# Set memory usage
if [[ -z ${MEMORY_TO_USE} ]]; then
    echo "No memory limit set by user. Defaulting to 16GB"
    MEMORY_TO_USE="16"
fi
###################################################################################################

###################################################################################################
# ENSURE THAT ALL REQUIRED SOFTWARE IS INSTALLED
###################################################################################################
command -v blastn > /dev/null || \
{ echo -e "This program requires 'blastn'. \n" \
        "Please install from NCBI website or NCBI github and retry. \n" \
        "Exiting..."; exit 2
    }
###################################################################################################

###################################################################################################
# PROCESS INPUT NAME SO IT'S EASIER TO HANDLE
###################################################################################################
# Create names for BLAST output file: first, truncate file path, leaving just the filename itself
SEQDUMP_FILE=${SEQDUMP##*/}

# Next, eliminate file extension, giving a cleaner name for blast
BLAST_NAME_SEQDUMP=${SEQDUMP_FILE%.*}
###################################################################################################

###################################################################################################
# READ ALL INPUTS BACK TO USER
###################################################################################################
# Create log file
LOG_FILE=./${BLAST_NAME_SEQDUMP}.cleanup_blast.log
touch ${LOG_FILE}

# Copy initial launch command into the log
echo -e "\ncleanup_blast was launched with the following command: \n    '${BASH_SOURCE} ${ALL_PARAMETERS}' \n" \
        "at: `date`" | tee -a ${LOG_FILE}

# Read inputs back to the user and store them in the log
echo -e "\n" \
        "Input fasta (seqdump) file provided: ${SEQDUMP} \n" \
        "Database used: ${PATH_TO_NT_DB} \n" \
        "e-value: ${E_VALUE} \n" \
        "Blast program: nucleotide-blast > ${BLAST_TASK} \n" \
        "Number of processors to use: ${NUM_THREADS} \n" \
        "Memory limit: ${MEMORY_TO_USE}GB \n"| tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# Run nucleotide blast
###################################################################################################
blastn \
-task ${BLAST_TASK} \
-db ${PATH_TO_NT_DB} \
-query ${SEQDUMP} \
-out ./${BLAST_NAME_SEQDUMP}.cleanup.results.txt \
-evalue ${E_VALUE} \
-num_threads ${NUM_THREADS} \
-outfmt "6 qseqid evalue stitle" \
-max_target_seqs 10000000 \
-max_hsps 1

# max_hsps only restricts to best level match PER VIRUS;
# if one read matches to multiple viruses, there will be multiple viruses reported

# Sort the reads, first by name then by evalue with lowest (strongest) on top; then, deduplicate
sort -k1,1 -k2g,2g ${BLAST_NAME_SEQDUMP}.cleanup.results.txt  | sort -k1,1 -u > \
    ${BLAST_NAME_SEQDUMP}.cleanup.results.top-hits.txt

# Overwrite the original hits file with the deduplicated results files
mv ./${BLAST_NAME_SEQDUMP}.cleanup.results.top-hits.txt ./${BLAST_NAME_SEQDUMP}.cleanup.results.txt
###################################################################################################

###################################################################################################
# OUTPUT LOGS
###################################################################################################
# Make a token to indicate the job finished correctly
echo -e "Finished nucleotide BLAST (${BLAST_TASK}), using file '${BLAST_NAME_SEQDUMP}' to query against \n" \
        "blast database '${PATH_TO_NT_DB}' at: \n`date` \n" | tee -a ${LOG_FILE}

# Print number of sequences in the input file:
echo -e "Number of sequences in original input file: "\
        "`grep -c "^>" ${SEQDUMP}`" | tee -a ${LOG_FILE}

# Print number of hits
echo -e "Number of hits in cleaned output hits list: " \
        "`wc -l ./${BLAST_NAME_SEQDUMP}.cleanup.results.txt | awk '{$1=$1};1' | cut -d " " -f 1` \n" | \
        tee -a ${LOG_FILE}

###################################################################################################

###################################################################################################
# FIND READS THAT MAPPED INITIALLY, BUT FELL OUT DURING THIS CLEANUP STEP
###################################################################################################
# Create temp files with reads in input vs. reads in output
INPUT_READS="input_reads.${BLAST_NAME_SEQDUMP}"
OUTPUT_READS="output_reads.${BLAST_NAME_SEQDUMP}"

grep "^>" ${SEQDUMP} | sed 's/>//g' | cut -d " " -f 1 | sort > ${INPUT_READS}
cut -f 1 ${BLAST_NAME_SEQDUMP}.cleanup.results.txt | cut -d " " -f 1 | sort > ${OUTPUT_READS}

# Find the reads in one but not the other; extract those reads from the input
seqtk subseq ${SEQDUMP} <(comm -3 ${INPUT_READS} ${OUTPUT_READS}) > \
    ${BLAST_NAME_SEQDUMP}.cleanup.dropped-reads.fasta

# Remove the temporary files
rm ${INPUT_READS}
rm ${OUTPUT_READS}

# Print number of dropped out reads
echo -e "Number of reads that dropped out of analysis during this cleaning step: " \
        "`grep -c "^>" ${BLAST_NAME_SEQDUMP}.cleanup.dropped-reads.fasta | \
          cut -d " " -f 1` \n" | \
        tee -a ${LOG_FILE}

# Save the names of those dropped out reads to the log
echo -e "Names of any reads that dropped out during the analysis: \n" >> ${LOG_FILE}

grep "^>" ${BLAST_NAME_SEQDUMP}.cleanup.dropped-reads.fasta | \
    tr -d ">" >> \
    ${LOG_FILE}
###################################################################################################

