FROM --platform=linux/amd64 ubuntu:20.04

ENV PYTHONPATH=/home/fpt/bin

WORKDIR /home/fpt

COPY . .

RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get install -y libgmp-dev \
    && apt-get install -y libmpfr-dev \
    && apt-get install -y libmpc-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

RUN . /opt/conda/etc/profile.d/conda.sh && \
    conda create -n "fpt_env" python=3.11.5 -y && \
    conda activate fpt_env && \
    conda install pip -y && \
    pip install -r minimal_requirments.txt && \
    conda init bash && \
    conda activate fpt_env && \ 
    make dependencies_docker && \
    make ladders_docker

ENV PATH=$CONDA_DIR/bin:$PATH

EXPOSE 8888


