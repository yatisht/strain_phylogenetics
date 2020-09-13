FROM ubuntu:18.04
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
USER root
WORKDIR HOME
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git \
    ca-certificates \
    sudo 
RUN git clone https://github.com/yatisht/strain_phylogenetics.git 
WORKDIR strain_phylogenetics
RUN ./installUbuntu.sh 
