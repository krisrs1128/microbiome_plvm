FROM r-base:3.3.3

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('plyr')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('feather')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('tidyr')"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "devtools::install_github('krisrs1128/ggscaffold')"

RUN apt-get update
RUN apt-get install python3-pip
RUN pip3 install luigi