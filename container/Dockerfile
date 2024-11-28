FROM rocker/tidyverse
RUN apt-get update && apt-get install -y parallel
RUN install2.r --deps TRUE --error --skipinstalled --ncpus -1 \
    survival \
    JMbayes2 \
    nlme