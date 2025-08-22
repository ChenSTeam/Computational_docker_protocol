# Base image
FROM ubuntu:25.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

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

RUN conda install conda-forge::openmm openmmforcefields mdtraj pdbfixer matplotlib -y && \
    conda install conda-forge::libsqlite --force-reinstall -y && \
    pip install torch torchvision --index-url https://download.pytorch.org/whl/cu126

# Set working directory
WORKDIR /app
COPY ./scripts /app

# Default command = run Python
ENTRYPOINT ["conda", "run", "-n", "md", "python"]