#!/usr/bin

# conda activate MEASubsampling


#snakemake -s Snakefile -j 1 --latency-wait 60

snakemake -j 4 --latency-wait 60 --rerun-incomplete

# End of file
