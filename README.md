# Latent Variable Modeling for the Microbiome

This repository includes all code for simulations, data analysis, and figures in
"Latent Variable Modeling for the Microbiome" by Sankaran and Holmes. The high
level commands, corresponding to different sections, are given in
`src/run_all.sh`. To facilitate direct experimentation, we have also provided a
docker image (built using the Dockerfile, and available from
the [Stanford Digital Repository](http://purl.stanford.edu/mf740kp2295) with all
software requirements preinstalled and this repo copied into the image's
`/home/` directory. The image can be loaded using `docker load --input
microbiome_plvm.tar`, and analysis can be reproduced by opening a container
`docker run -it microbiome_plvm bash`, navigating to `/home/src/` and executing
lines from `run_all.sh`.
