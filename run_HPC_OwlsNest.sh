#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > MESSI_dag.pdf

mkdir -p logs

snakemake \
      -s Snakefile \
      --cluster-config cluster.json \
      --cluster "qsub -V -N MESSI -l nodes={cluster.nodes}:ppn={cluster.ppn} -l walltime=48:00:00 -q {cluster.name} -e logs -o logs" \
      --jobs 8 all \
      --keep-going \
      --reason \
      --rerun-incomplete \
      --latency-wait 120 

# End of file 
