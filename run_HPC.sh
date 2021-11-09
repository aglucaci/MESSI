#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

snakemake \
      -s Snakefile \
      --cluster-config cluster.yml \
      --cluster "qsub -V -N MESSI -l nodes={cluster.nodes}:ppn={cluster.ppn} -l walltime=120:00:00 -q {cluster.name} -e logs -o logs" \
      --jobs 1 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 \
      --use-conda

# End of file 
