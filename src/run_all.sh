source .env

## Runs simulation and generates figures for the NMF experiment
#Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_expers.R
#Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_vis.R

cd $MICROBIOME_PLVM_DIR/src/sim/lda/pipeline
pwd
Rscript ../Rscript/build_config.R
python3 pipeline.py LDAExperiment --local-scheduler --workers=8
cd $MICROBIOME_PLVM_DIR
