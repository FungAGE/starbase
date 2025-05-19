# Select base image
FROM condaforge/miniforge3
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

# System dependencies and redis installation
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables wget redis-server && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy environment file and create conda environment 
COPY environment.yaml .
RUN conda env create -f environment.yaml && \
    conda clean -afy

# Set conda environment to activate by default
ENV PATH=/opt/conda/envs/starbase/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase

# Install Node.js, npm, and blasterjs
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs

# Create cache directory
RUN mkdir -p $HOME/src/database/db/cache && \
    chmod -R 777 $HOME/src/database/db/cache && \
    chown -R $USER:$USER $HOME

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