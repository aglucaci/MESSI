# Molecular Evolution SubSampling Investigator (MESSI)

## Install conda environment 
> conda env create -f environment.yml

This will create a virtual environment called "MESSI" with the necessary dependencies.

You can also start this from the "Initialize.sh" file with ```bash Initialize.sh```

## Create tmux session

A helper script is included to create a tmux session titled "MESSI", you can enable this with ```bash MakeTmuxSession.sh```

## Configure your HPC environment

A yaml configuration file is included ```cluster.yml``` this species the name and ppn of your HPC environments job task queue to use.

## Configure your test dataset

A yaml configuration file is included ```config.yml``` this specifies the file names of your input alignment (required) and newick tree (required). These files should be placed in the ```Data/``` folder directly. Current design also requires that you specify your working directory "WD" and the output directory, typicall this is the "Results/" folder.

## Default Data

> Data/13-Datasets_Split/benchmark/
Contains the 13 empirical datasets (in NEXUS format) split out to their corresponding FASTA and Newick file. With accompanying scripts.

Place your alignment and newick tree under investigation directly into the "Data/" folder.
