## Why
This script was made to easily split a folder of fastqs into smaller parts and run the gms-arti pipeline on each chunk. 
The background to this was that the NGPc had problems tracking a very large number of jobs, as produced by hundreds of samples.#

## Known Issues
This will currently only work for the Illumina based workflow. The nanopore workflow requires a different input in terms of data structure, and it will not be handled 
by this script.
