FROM nvidia/cuda:12.1.0-base-ubuntu22.04

### Dockerfile for running basecalling from AWS Batch

RUN apt-get update && apt-get install -y --no-install-recommends \
  wget curl git ca-certificates build-essential nvidia-cuda-toolkit \
  libhdf5-dev libssl-dev libzstd-dev cmake autoconf automake samtools

WORKDIR /tmp

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
RUN chmod +x ./Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
RUN ./Miniconda3-py39_23.1.0-1-Linux-x86_64.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:$PATH
RUN conda init

RUN pip install duplex_tools
RUN pip install awscli

RUN useradd -ms /bin/bash ubuntu
USER ubuntu

WORKDIR /home/ubuntu
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.4.1-linux-x64.tar.gz
RUN tar -xvf dorado-0.4.1-linux-x64.tar.gz

RUN ./dorado-0.4.1-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0
RUN ./dorado-0.4.1-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mC@v2
RUN ./dorado-0.4.1-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0_6mA@v3

COPY --chown=ubuntu:ubuntu ./run_dorado.sh ./run_dorado.sh
RUN chmod +x /home/ubuntu/run_dorado.sh

RUN mkdir -p /home/ubuntu/basecalling
WORKDIR /home/ubuntu/basecalling

ENTRYPOINT ["/home/ubuntu/run_dorado.sh"]
