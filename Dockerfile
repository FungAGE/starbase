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

# System dependencies (as root)
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables wget redis-server build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Node.js, npm, and blasterjs (as root)
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs

# Switch to user for conda/mamba installation
USER $USER

# Install Mambaforge as the user (so files are owned by user from the start)
RUN wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $HOME/miniconda && \
    rm miniforge.sh && \
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate starbase" >> $HOME/.bashrc

ENV PATH=$HOME/miniconda/bin:$PATH

# Copy environment file and create conda environment using mamba
COPY --chown=$USER:$USER environment.yaml .
RUN mamba env create -y -f environment.yaml && \
    mamba clean -afy

# Set conda environment to activate by default
ENV PATH=$HOME/miniconda/envs/starbase/bin:$HOME/miniconda/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase

# Create required directories (as user, before copying app code)
RUN mkdir -p $HOME/src/database/db \
             $HOME/src/database/logs \
             $HOME/src/database/cache && \
    chmod -R 755 $HOME/src/database/logs \
                 $HOME/src/database/cache && \
    chmod -R 777 $HOME/src/database/db

# Add healthcheck
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/api/cache/status || exit 1

# Copy application code (changes most frequently, so do this last)
COPY --chown=$USER:$USER ./ ./
RUN chmod +x start-script.sh && \
    chmod +x start_celery.py && \
    chmod +x manage_celery.py

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]