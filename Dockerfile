FROM continuumio/miniconda:4.5.4

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY envs/linux_env.yml /environment.yml
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/pipeliner/bin:$PATH