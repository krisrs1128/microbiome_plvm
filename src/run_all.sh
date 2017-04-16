source .env

## Runs simulation and generates figures for the NMF experiment
Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_expers.R
Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_vis.R

## Run simulations and generate figures for LDA experiment
cd $MICROBIOME_PLVM_DIR/src/sim/lda/pipeline
Rscript ../Rscript/build_config.R
python3 pipeline.py LDAExperiment --local-scheduler --workers=2
Rscript ../lda_vis.R

cd $MICROBIOME_PLVM_DIR
