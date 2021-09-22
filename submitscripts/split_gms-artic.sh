#!/bin/bash

### #CLI
#Defaults
CHUNKSIZE=50

#Usage instructions on no args
if  [ "$#" == 0 ]; then
    SCRIPTNAME=`basename "$0"`
    echo >&2 -e "Usage: $SCRIPTNAME 
    \t-f /path/to/fastq/dir
    \t-t /path/to/temp/dir
    \t-a /path/to/gms-artic
    \t-o /path/to/output/dir
    \t-p Prefix for output files
    \t-c Chunk size to use (Default: $CHUNKSIZE)" 
    exit 0
fi

#Read in all flags
while getopts ":f:t:a:o:p:c:" opt; do
    case $opt in
	f)
	    FASTQDIR=$OPTARG
	    ;;
	t)
	    TEMPFASTQ=$OPTARG
	    ;;
	a)
	    ARTICDIR=$OPTARG
	    ;;
	o)
	    OUTPUTDIR=$OPTARG
	    ;;
	p)
	    PREFIX=$OPTARG
	    ;;
	c)
	    CHUNKSIZE=$OPTARG
	    ;;
	\?)
	    echo >&2 "Invalid option: -$OPTARG" >&2
	    exit 1
	    ;;
	:)
	    echo >&2 "Option -$OPTARG requires an argument." >&2
	    exit 1
	    ;;
    esac
done

### Sanity checks
# Does tempdir exist? 
if [ ! -d $TEMPFASTQ ]; then
    echo "Directory $TEMPFASTQ DOES NOT exists. Please create it."
    exit 1
fi

# Is temp dir empty?
if [ "$(ls -A $TEMPFASTQ)" ]; then
    echo "$TEMPFASTQ is not empty, please correct this."
    exit 1
fi

# Are there any FASTQ files in the FASTQDIR?
if [ $(ls ${FASTQDIR}/*.fastq.gz | wc -l) -lt 1 ]; then
    echo "No files with ending .fastq.gz found in $FASTQDIR"
    exit 1
fi

#Is BASH version 4 or higher (To use negative array indecies
if [ "${BASH_VERSINFO:-0}" -lt 4 ]; then
    echo "BASH version < 4, won't work!"
    exit 1
fi

#Is there a main.nf in ARTICDIR
if [ ! -f "${ARTICDIR}/main.nf" ]; then
    echo "Could not find a main.nf in $ARTICDIR."
    exit 1
fi

#Set up some counters
COUNTER=0
RESCOUNT=1

#Load in all files to process into an array
FASTQS=($FASTQDIR/*fastq.gz)

for FASTQ in ${FASTQS[@]}; do
    let COUNTER=COUNTER+1
    rsync -P $FASTQ ${TEMPFASTQ}/ 

    #Run CHUNKSIZE num of samples at once
    if [[ $COUNTER -gt $CHUNKSIZE-1 ]]; then
	#Run pipeline
	nextflow run ${ARTICDIR}/main.nf -profile singularity,sge --illumina --prefix $PREFIX --directory $TEMPFASTQ --outdir ${OUTPUTDIR}-$RESCOUNT

	#Reset counters
	let RESCOUNT=RESCOUNT+1
	COUNTER=0
	
	#Remove temp data
	rm $TEMPFASTQ/*.fastq.gz
    fi

    #Run again for all remaining samples
    if [ $FASTQ == ${FASTQS[-1]} ] && [ $COUNTER -gt 0 ]; then
	nextflow run ${ARTICDIR}/main.nf -profile singularity,sge --illumina --prefix $PREFIX --directory $TEMPFASTQ --outdir ${OUTPUTDIR}-$RESCOUNT
	rm $TEMPFASTQ/*.fastq.gz
    fi


done
