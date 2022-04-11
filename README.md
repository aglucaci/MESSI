# Molecular Evolution SubSampling Investigator (MESSI)

## Install conda environment 
> conda env create -f environment.yml

This will create a virtual environment called "MESSI" with the necessary dependencies.

You can also start this from the "Initialize.sh" file with ```bash Initialize.sh```

## Create tmux session

A helper script is included to create a tmux session titled "MESSI", you can enable this with ```bash MakeTmuxSession.sh```

## Configure your HPC environment

A JSON configuration file is included ```cluster.json``` this species the name and ppn of your HPC environments job task queue to use.

## Configure your test dataset

A YAML configuration file is included ```config.yml``` this specifies the file names of your input alignment (required) and newick tree (required). You should also specify your dataset label via the "Tag" parameter, this will be the subdirectory of the "Results" folder where your output will be stored. These files should be placed in the ```Data/``` folder directly. Current design also requires that you specify your working directory "WD" and the output directory, typicall this is the "Results/" folder.

## Default empirical data

> Data/13-Datasets_Split/benchmark/

Contains the 13 empirical datasets (in NEXUS format) split out to their corresponding FASTA and Newick file. With accompanying scripts ```seqmagick_convert.sh``` which is used to convert NEXUS formatted files to FASTA. Along with a custom script ```nexus_extract_tree.py``` used to extract the Newick tree from a nexus file.

Note: Place your alignment and newick tree under investigation directly into the "Data/" folder.

## Simulated data

The folder ```Data/Simulations/``` holds preliminary simulated data. It also contains a bash script ```CreateSims.sh``` used to create our simulated data.
