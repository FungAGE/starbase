# Select base image
FROM python:3.9
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="STARBASE is a database and toolkit for exploring large transposable elements in Fungi"

ARG IPSTACK_API_KEY
ARG MAINTENANCE_TOKEN

# Create variables for user name, home directory, and secrets
ENV USER=starbase
ENV HOME=/home/$USER
ENV IPSTACK_API_KEY=$IPSTACK_API_KEY
ENV MAINTENANCE_TOKEN=$MAINTENANCE_TOKEN

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory
WORKDIR $HOME/

# System dependencies
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables wget redis-server build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Mambaforge
RUN wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $HOME/miniconda && \
    rm miniforge.sh && \
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate starbase" >> $HOME/.bashrc

ENV PATH=$HOME/miniconda/bin:$PATH

# Copy environment file and create conda environment using mamba
COPY environment.yaml .
RUN mamba env create -f environment.yaml && \
    mamba clean -afy

# Set conda environment to activate by default
ENV PATH=$HOME/miniconda/envs/starbase/bin:$HOME/miniconda/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase

# Install Node.js, npm, and blasterjs
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs

# Create db directory
RUN mkdir -p $HOME/src/database/db && \
    chown -R $USER:$USER $HOME && \
    chmod -R 777 $HOME/src/database/db

# Add healthcheck
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/api/cache/status || exit 1

# Copy application code (changes most frequently, so do this last)
COPY ./ ./
RUN chmod +x start-script.sh && \
    # Ensure all directories and files are owned by starbase user
    chown -R $USER:$USER $HOME/src

# Switch to user
USER $USER

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]