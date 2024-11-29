FROM rocker/tidyverse
LABEL "org.opencontainers.image.url"="https://github.com/anand-imcm/joint-model-analysis"
LABEL "org.opencontainers.image.title"="joint-model-analysis"
LABEL "org.opencontainers.image.description"="Tools for Joint Modelling analysis"
LABEL "org.opencontainers.image.source"="https://github.com/anand-imcm/joint-model-analysis"
RUN apt-get update && apt-get install -y parallel
COPY ./scripts /scripts
RUN install2.r --deps TRUE --error --skipinstalled --ncpus -1 \
    survival \
    JMbayes2 \
    nlme