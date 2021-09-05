#!/usr/bin

# conda activate MEASubsampling


snakemake -s Snakefile -j 1 --latency-wait 60
