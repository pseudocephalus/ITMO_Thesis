# R environment for scripts/SCoRe.R and scripts/GWAS.R
FROM rocker/tidyverse:4.4.2

# Ancestry clustering package (GitHub only)
RUN R -e 'remotes::install_github("alexloboda/SVDFunctions")'

WORKDIR /work