FROM continuumio/miniconda
COPY envs/linux_env.yml /
RUN conda env create -f /linux_env.yml && conda clean -a
ENV PATH /opt/conda/envs/pipeliner/bin:$PATH