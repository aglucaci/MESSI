#!/bin/bash
#PBS -N SimMG94

#@Usage: bash CreateSims.sh

# Hyphy related
#HYPHYMPI="/home/aglucaci/hyphy-develop/HYPHYMPI"
#HYPHY="/home/aglucaci/hyphy-develop/HYPHYMP"
#HYPHY=hyphy
#RES="/home/aglucaci/hyphy-develop/res"
#SIM_MG94="/home/aglucaci/hyphy-analyses/SimulateMG94/SimulateMG94.bf"

HYPHY=hyphy
SIM_MG94=/home/tuk13147/work/hyphy-analyses/SimulateMG94/SimulateMG94.bf

# Simulations related inputs
#TREE="/home/aglucaci/SubsamplingSequences/data/TP53_refseq_cds_CodonAligned.nwk"
#TREE=/home/tuk13147/work/MESSI/Data/TP53_refseq_cds_CodonAligned.nwk
TREE=/home/tuk13147/work/MESSI/Data/HIV_RT.fasta.treefile

# Output
#OUTPUT="/home/aglucaci/SubsamplingSequences/data/Sims"
OUTPUT=/home/tuk13147/work/MESSI/Data/Simulations
#OUTPUT=/home/tuk13147/work/MESSI/Data
mkdir -p $OUTPUT

#NP=16

# #####################################################################
# Main
# #####################################################################

## Test phase
#hyphy SimulateMG94.bf --tree CD2.nwk --omega 0.25 --sites 300 --replicates 10 --output data/example1

## First set of sims
#$HYPHY LIBPATH=$RES $SIM_MG94 --tree $TREE --omega 1 --sites 300 --replicates 10 --output $OUTPUT"/p53_simulations_o1_s300_r10"
#$HYPHY LIBPATH=$RES $SIM_MG94 --tree $TREE --omega 0.1 --sites 300 --replicates 10 --output $OUTPUT"/p53_simulations_o01_s300_r10"
#$HYPHY LIBPATH=$RES $SIM_MG94 --tree $TREE --omega 10 --sites 300 --replicates 10 --output $OUTPUT"/p53_simulations_o10_s300_r10"

## Second set of sims
# gamma distributed srv
#mpirun -np $NP $HYPHYMPI LIBPATH=$RES $SIM_MG94 --tree $TREE --omega 0.1 --sites 300 --replicates 100 --site-variation gamma --output $OUTPUT"/TP53_simulations_o01_s300_r100_srv_gamma"
#mpirun -np $NP $HYPHYMPI LIBPATH=$RES $SIM_MG94 --tree $TREE --omega 1 --sites 300 --replicates 100 --site-variation gamma --output $OUTPUT"/TP53_simulations_o1_s300_r100_srv_gamma"
#mpirun -np $NP $HYPHYMPI LIBPATH=$RES $SIM_MG94 --tree $TREE --omega 10 --sites 300 --replicates 100 --site-variation gamma --output $OUTPUT"/TP53_simulations_o10_s300_r100_srv_gamma"

# Last run -- TP53
#$HYPHY $SIM_MG94 --tree $TREE --omega 1 --sites 500 --replicates 1 --output $OUTPUT"/TP53_simulations_o1_s500_r1"
#$HYPHY $SIM_MG94 --tree $TREE --omega 0.1 --sites 500 --replicates 1 --output $OUTPUT"/TP53_simulations_o01_s500_r1"
#$HYPHY $SIM_MG94 --tree $TREE --omega 10 --sites 500 --replicates 1 --output $OUTPUT"/TP53_simulations_o10_s500_r1"

# HIV_RT
$HYPHY $SIM_MG94 --tree $TREE --omega 1 --sites 300 --replicates 1 --output $OUTPUT"/HIV-RT_simulations_o1_s300_r1"
$HYPHY $SIM_MG94 --tree $TREE --omega 0.1 --sites 300 --replicates 1 --output $OUTPUT"/HIV-RT_simulations_o01_s300_r1"
$HYPHY $SIM_MG94 --tree $TREE --omega 10 --sites 300 --replicates 1 --output $OUTPUT"/HIV-RT_simulations_o10_s300_r1"

# Split from FNA to pasta
yourfilenames=`ls $OUTPUT/*.replicate.*`

for f in $yourfilenames; do
    echo "Splitting fasta and newick: "$f
    if [ $f == "$f".nwk ]; then
        continue
    fi

    if [ $f == "$f".fasta ]; then
        continue
    fi

    tail -n 1 $f > "$f".newick
    head -n-1 $f > "$f".fasta
done



exit 0

# #####################################################################
# End of file
# #####################################################################
