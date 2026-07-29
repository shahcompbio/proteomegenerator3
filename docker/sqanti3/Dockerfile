FROM quay.io/biocontainers/sqanti3:6.0.1--hdfd78af_0

USER root

# Install RColorConesa (missing from biocontainers image, needed for filter report)
RUN R -e "install.packages('RColorConesa', repos='https://cran.rstudio.com/')"

# Set timezone to avoid R issues
ENV TZ=UTC

USER bioc
