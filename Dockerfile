FROM condaforge/miniforge3:latest
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
RUN useradd -m -u 1001 $USER

# Set working directory
WORKDIR $HOME/

# System dependencies and supercronic installation
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables wget redis-server && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    # Install supercronic
    curl -fsSLO "$SUPERCRONIC_URL" && \
    echo "${SUPERCRONIC_SHA1SUM}  ${SUPERCRONIC}" | sha1sum -c - && \
    chmod +x "$SUPERCRONIC" && \
    mv "$SUPERCRONIC" "/usr/local/bin/supercronic"

# Copy environment file and create conda environment
COPY environment.yaml .
RUN conda env create -f environment.yaml && \
    conda clean -afy

# Copy application code (changes most frequently, so do this last)
COPY ./ ./

# Create necessary directories and set permissions
RUN mkdir -p $HOME/src/database/db/cache /var/run/crond /var/log/cron && \
    chmod +x start-script.sh && \
    # Ensure all directories and files are owned by starbase user
    chown -R $USER:$USER $HOME && \
    chmod -R 755 $HOME && \
    chmod -R 777 $HOME/src/database/db/cache && \
    chmod -R 777 /var/run/crond /var/log/cron

# Set conda environment to activate by default
SHELL ["conda", "run", "-n", "starbase", "/bin/bash", "-c"]

# Install Node.js, npm, and blasterjs
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs

# Create cron and cache directories and set up logging
RUN mkdir -p /var/run/crond /var/log/cron $HOME/cron $HOME/src/database/db/cache && \
    touch /var/log/cron/cron.log && \
    # Create crontab file for supercronic
    echo "0 * * * * cd $HOME && python -m src.utils.telemetry update_ip_locations >> $HOME/cron/cron.log 2>&1" > $HOME/cron/crontab && \
    echo "*/15 * * * * cd $HOME && curl -X POST http://localhost:8000/api/refresh-telemetry >> $HOME/cron/cron.log 2>&1" >> $HOME/cron/crontab && \
    # Add cache check to crontab
    echo "*/5 * * * * curl -f http://localhost:8000/api/cache/status || curl -X POST http://localhost:8000/api/cache/refresh" >> $HOME/cron/crontab && \
    # Add Celery beat check (hourly check to make sure it's running)
    echo "0 * * * * if ! pgrep -f 'celery -A src.config.celery_config:celery beat' > /dev/null; then cd $HOME && restart_celery_beat >> $HOME/cron/cron.log 2>&1; fi" >> $HOME/cron/crontab && \
    # Set permissions for all directories
    chown -R $USER:$USER $HOME /var/run/crond /var/log/cron && \
    chmod -R 777 $HOME/src/database/db/cache

# Add healthcheck
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/api/cache/status || exit 1

# Switch to user
USER $USER

# Always activate conda environment when container starts
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "starbase"]
CMD ["./start-script.sh"]