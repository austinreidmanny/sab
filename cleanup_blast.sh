#!/bin/bash

###################################################################################################
# ENSURE THE SCRIPT IS CALLED CORRECTLY
###################################################################################################
usage() { echo -e "The objective of this script is to clean up the results of stand-alone-blast (sab). \n\n" \
                  "This program takes in a FASTA (given as output by sab) of reads that mapped " \
                  "to the previous reference and searches those again, this time against the \n" \
                  "entire NCBI nonredundant (nr) database. The output of this script is a \n" \
                  "tab-delimited text file with the following fields: \n" \
                  "(1) name of the read (2) taxon of best hit in NCBI (3) e-value score of match \n\n" \
                  "Usage: $0 -i input.fasta [options] \n\n" \
                  "Output: input.cleanup.results.txt \n\n" \
                  "Optional parameters: \n" \
                        "-e (evalue, e.g. 100, 1, or 1e-99; [default = 1e-9]) \n" \
                        "-m (maximum amount of memory to use [in GB]; [default=16] ) \n" \
                        "-p (path to directory for NCBI nt database; [default='~/Documents/Research/sra/blastdbs/'] ) \n" \
                        "-n (sets nucleotide program to blastn; [default= dc-megablast] ) \n" \
                        "-g (sets nucleotide program to megablast; [default= dc-megablast] ) \n\n" \
                      "Example of a complex run: \n" \
                      "$0 -i input.fasta -e 1e-3 -m 26 -g \n\n" \
                      "Exiting program. Please retry with corrected parameters..." >&2 && exit 1
        }
###################################################################################################

###################################################################################################
# TAKE IN THE USER-PROVIDED PARAMETERS
###################################################################################################
while getopts "i:e:m:p:n*:" arg; do
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
    echo -e "ERROR: No input detected. \n\n"
    usage
fi

# If e-value wasn't provided by user, then set it to 1e-9
if [[ -z ${E_VALUE} ]]; then
    E_VALUE="1e-9"
fi

# If -n (blastn) or -g (megablast) flags were not given by user, default to dc-megablast
if [[ -z "${BLAST_TASK}" ]]; then
    BLAST_TASK="dc-megablast"
fi

# If path to NCBI nt database was not given, give default path
if [[ -z "${PATH_TO_NT_DB}" ]]; then
    PATH_TO_NT_DB="${HOME}/Documents/Research/sra/blastdbs/nt"
fi
###################################################################################################

###################################################################################################
# Set up number of CPUs to use (use all available) and RAM (can be set by user, defaults to 16GB)
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
LOG_FILE=./${BLAST_NAME_SEQDUMP}.cleanup.log
touch ${LOG_FILE}

# Copy initial launch command into the log
echo -e "sab was launched with the following command: \n $0 $@ \n" > ${LOG_FILE}

# Read inputs back to the user and store them in the log
echo -e "\n" \
        "Input fasta (seqdump) file provided: ${SEQDUMP} \n" \
        "Database used: NCBI NT databse (located at: ${PATH_TO_NT_DB})"
        "e-value: ${E_VALUE} \n" \
        "Blast program: nucleotide-blast > ${BLAST_TASK} \n" \
        "Number of processors to use: ${NUM_THREADS} \n" \
        "Memory limit: ${MEMORY_TO_USE}GB \n\n"| tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# Run nucleotide blast
###################################################################################################
blastn \
-task ${BLAST_TASK} \
-db ${PATH_TO_NT_DB}\
-query ${SEQDUMP} \
-out ./${BLAST_NAME_SEQDUMP}.cleanup.results.txt \
-evalue ${E_VALUE} \
-num_threads ${NUM_THREADS} \
-outfmt "6 qseqid sseqid evalue" \
-max_target_seqs 10000000
###################################################################################################

###################################################################################################
# OUTPUT LOGS
###################################################################################################
# Make a token to indicate the job finished correctly
echo "Finished nucleotide BLAST (${BLAST_TASK}), using ${SEQDUMP} to query against NCBI NT database" | \
tee ${LOG_FILE}

# Indicate time of completion
echo "Job finished at" && date | tee ${LOG_FILE}
###################################################################################################
