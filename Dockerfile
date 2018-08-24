FROM continuumio/miniconda2:latest

WORKDIR /

# Install requirements
COPY envs/linux_env.yml /environment.yml
RUN conda config --add channels conda-forge \
    && conda env create -n pipeliner -f environment.yml \
    && rm -rf /opt/conda/pkgs/*

# Install
RUN chown -R 777 /*

# Activate environment
ENV PATH /opt/conda/envs/pipeliner/bin:$PATH