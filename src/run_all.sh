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
cd $MICROBIOME_PLVM_DIR/src/sim/lda/pipeline
Rscript ../Rscript/build_config.R
python pipeline.py LDAExperiment --local-scheduler --workers=2
Rscript ../lda_vis.R

## Run simulations and generate figures for Unigram experiment
cd $MICROBIOME_PLVM_DIR/src/sim/unigram/pipeline
Rscript ../Rscript/build_config.R
python pipeline.py UnigramExperiment --local-scheduler --workers=2
# Rscript ../unigram_vis.R

cd $MICROBIOME_PLVM_DIR

## Run and save for antibiotics application to LDA experiment
cd $MICROBIOME_PLVM_DIR/src/antibiotics-study
chmod +x lda.R unigram.R posterior_checks.R
./lda.R
./unigram.R
./posterior_checks.R
./lda.R --subject 'D'
./unigram.R --subject 'D'
./posterior_checks.R --subject 'D'
./lda.R --subject 'E'
./unigram.R --subject 'E'
./posterior_checks.R --subject 'E'

cd $MICROBIOME_PLVM_DIR
