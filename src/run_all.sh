source .env

## Runs simulation and generates figures for the NMF experiment
Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_expers.R
Rscript $MICROBIOME_PLVM_DIR/src/sim/nmf/nmf_vis.R
