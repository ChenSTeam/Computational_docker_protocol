FROM ubuntu:25.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda.sh && \
    bash ./miniconda.sh -b -p $CONDA_DIR && \
    rm ./miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# Create conda environment from environment.yml
RUN conda tos accept && \
    conda create -n md python=3.12 -y

# Activate environment by default
SHELL ["conda", "run", "-n", "md", "/bin/bash", "-c"]

RUN conda install conda-forge::openmm openmmforcefields mdtraj pdbfixer matplotlib gputil cudatoolkit cudnn -y && \
    conda install conda-forge::libsqlite conda-forge::seaborn --force-reinstall -y

# Set working directory
WORKDIR /app
COPY ./scripts /app

# Default command = run Python
ENTRYPOINT ["conda", "run", "-n", "md", "python", "md_rama_protocol.py"]