#!/bin/sh
# Central location to run all simulations and data analysis appearing in "Latent
# Variable Modeling for the Microbiome"
#
# author: kriss1@stanford.edu

# define $MICROBIOME_PLVM_DIR environmental variable
source .env

## Runs simulation and generates figures for the NMF experiment
Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_expers.R
Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_vis.R

## Run simulations and generate figures for LDA experiment
# cd $MICROBIOME_PLVM_DIR/src/sim/lda/pipeline
Rscript ../Rscript/build_config.R
python3 pipeline.py LDAExperiment --local-scheduler --workers=2
Rscript ../lda_vis.R

cd $MICROBIOME_PLVM_DIR

## Run and save for antibiotics application to LDA experiment
cd $MICROBIOME_PLVM_DIR/src/antibiotics-study
Rscript lda.R
Rscript unigram.R

cd $MICROBIOME_PLVM_DIR
