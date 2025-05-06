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

# variables for supercronic installation (a cron for containers)
ENV SUPERCRONIC_URL=https://github.com/aptible/supercronic/releases/download/v0.2.24/supercronic-linux-amd64 \
    SUPERCRONIC=supercronic-linux-amd64 \
    SUPERCRONIC_SHA1SUM=6817299e04457e5d6ec4809c72ee13a43e95ba41

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory
WORKDIR $HOME/

# System dependencies, conda, and supercronic installation (combined to reduce layers)
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables wget redis-server && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    # Install supercronic
    curl -fsSLO "$SUPERCRONIC_URL" && \
    echo "${SUPERCRONIC_SHA1SUM}  ${SUPERCRONIC}" | sha1sum -c - && \
    chmod +x "$SUPERCRONIC" && \
    mv "$SUPERCRONIC" "/usr/local/bin/supercronic"

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p $HOME/miniconda && \
    rm miniconda.sh && \
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate starbase" >> $HOME/.bashrc

ENV PATH=$HOME/miniconda/bin:$PATH

# Copy environment file and create conda environment 
COPY environment.yaml .
RUN conda env create -f environment.yaml && \
    conda clean -afy

# Set conda environment to activate by default
ENV PATH=$HOME/miniconda/envs/starbase/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase

# Install Node.js, npm, and blasterjs
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs

# Create cron and cache directories and set up logging (combined directory creation)
RUN mkdir -p /var/run/crond /var/log/cron $HOME/cron $HOME/src/database/db/cache && \
    touch /var/log/cron/cron.log && \
    # Create crontab file for supercronic
    echo "0 * * * * cd $HOME && python -m src.utils.telemetry update_ip_locations >> $HOME/cron/cron.log 2>&1" > $HOME/cron/crontab && \
    echo "*/15 * * * * cd $HOME && curl -X POST http://localhost:8000/api/refresh-telemetry >> $HOME/cron/cron.log 2>&1" >> $HOME/cron/crontab && \
    # Add cache check to crontab
    echo "*/5 * * * * curl -f http://localhost:8000/api/cache/status || curl -X POST http://localhost:8000/api/cache/refresh" >> $HOME/cron/crontab && \
    # Set permissions for all directories
    chown -R $USER:$USER $HOME /var/run/crond /var/log/cron && \
    chmod -R 777 $HOME/src/database/db/cache

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