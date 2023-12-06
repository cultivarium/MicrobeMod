FROM python:3.9

RUN apt-get update

# Install conda
RUN mkdir -p ~/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
RUN bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
RUN rm -rf ~/miniconda3/miniconda.sh
RUN ~/miniconda3/bin/conda init bash

# Create directory for storing data
RUN mkdir /files

# Setup microbemod
COPY bin ./bin
COPY MicrobeMod ./MicrobeMod
COPY setup.py ./
RUN cd MicrobeMod && python download_db.py

# Install deps
RUN ~/miniconda3/bin/conda config --append channels bioconda \
  --append channels conda-forge \
  --append channels nanoporetech \
  --append channels anaconda
RUN ~/miniconda3/bin/conda install -c conda-forge libcxx
RUN ~/miniconda3/bin/conda install -c conda-forge libmamba
RUN ~/miniconda3/bin/conda install -c bioconda -y blast hmmer prodigal
RUN ~/miniconda3/bin/conda install -c bioconda -y cath-tools
RUN ~/miniconda3/bin/conda install -c bioconda -y meme
RUN ~/miniconda3/bin/conda install -c nanoporetech -y modkit
RUN ~/miniconda3/bin/conda install -c anaconda biopython pandas

COPY docker_microbemod.sh ./microbemod.sh

ENTRYPOINT ["bash", "./microbemod.sh"]
