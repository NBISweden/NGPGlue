#!/bin/bash

#CLI
FASTQDIR=$1
TEMPFASTQ=$2
ARTICDIR=$3
ARTICOUTPUT=$4
ARTICPREFIX=$5
CHUNKSIZE=${6-120}

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

#Set up some counters
COUNTER=0
RESCOUNT=1

#Load in al lfiles to process into an array
FASTQS=($FASTQDIR/*fastq.gz)

for FASTQ in ${FASTQS[@]}; do
    let COUNTER=COUNTER+1
    rsync -P $FASTQ ${TEMPFASTQ}/ 

    #Run CHUNKSIZE num of samples at once
    if [[ $COUNTER -gt $CHUNKSIZE-1 ]]; then
	#Run pipeline
	nextflow run ${ARTICDIR}/main.nf -profile singularity,sge --illumina --prefix $ARTICPREFIX --directory $TEMPFASTQ --outdir ${ARTICOUTPUT}-$RESCOUNT

	#Reset counters
	let RESCOUNT=RESCOUNT+1
	COUNTER=0
	
	#Remove temp data
	rm $TEMPFASTQ/*.fastq.gz
    fi

    #Run again for all remaining samples
    if [ $FASTQ == ${FASTQS[-1]} ] && [ $COUNTER -gt 0 ]; then
	nextflow run ${ARTICDIR}/main.nf -profile singularity,sge --illumina --prefix $ARTICPREFIX --directory $TEMPFASTQ --outdir ${ARTICOUTPUT}-$RESCOUNT
	rm $TEMPFASTQ/*.fastq.gz
    fi


done

