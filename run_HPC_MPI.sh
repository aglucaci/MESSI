#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

snakemake \
      -s Snakefile_MPI \
      --cluster-config cluster.yml \
      --cluster "qsub -V -N MESSI_MPI -l nodes={cluster.nodes}:ppn={cluster.ppn} -l walltime=120:00:00 -q {cluster.name} -e logs -o logs" \
      --jobs 8 all \
      --keep-going \
      --reason \
      --rerun-incomplete \
      --latency-wait 120 \
      --use-conda

# End of file 
