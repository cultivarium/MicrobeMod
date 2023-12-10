FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install wget python3 zlib1g-dev git build-essential -y

# Install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
RUN chmod +x ./Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
RUN ./Miniconda3-py39_23.1.0-1-Linux-x86_64.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:$PATH
RUN conda init

# Install deps
RUN wget https://meme-suite.org/meme/meme-software/5.5.5/meme-5.5.5.tar.gz
RUN tar -xvf meme-5.5.5.tar.gz
RUN cd meme-5.5.5 && ./configure --enable-build-libxml2 --enable-build-libxslt
RUN cd meme-5.5.5 && make
RUN mv meme-5.5.5/src/streme /usr/bin/
RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.2.3/modkit_v0.2.3_centos7_x86_64.tar.gz
RUN tar -xvf modkit_v0.2.3_centos7_x86_64.tar.gz
RUN mv dist/modkit /usr/bin/

RUN conda install -c conda-forge -c bioconda -y blast hmmer prodigal cath-tools

# Setup microbemod
RUN git clone https://github.com/cultivarium/MicrobeMod/
RUN cd MicrobeMod/MicrobeMod && python download_db.py
RUN cd MicrobeMod && pip install .


ENTRYPOINT ["MicrobeMod"]