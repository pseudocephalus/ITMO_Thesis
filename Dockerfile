# R environment for scripts/SCoRe.R and scripts/GWAS.R
FROM rocker/tidyverse:4.4.2

# Ancestry clustering package (GitHub only)
RUN R -e 'remotes::install_url("https://github.com/alexloboda/SVDFunctions/archive/86acf528b08f5ec7523a8d92aeba19e1f7dfe4d9.tar.gz")'

WORKDIR /work