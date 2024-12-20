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

# Install supercronic (a cron for containers)
ENV SUPERCRONIC_URL=https://github.com/aptible/supercronic/releases/download/v0.2.24/supercronic-linux-amd64 \
    SUPERCRONIC=supercronic-linux-amd64 \
    SUPERCRONIC_SHA1SUM=6817299e04457e5d6ec4809c72ee13a43e95ba41

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory
WORKDIR $HOME/

# System dependencies and supercronic installation (combined to reduce layers)
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables ncbi-blast+ hmmer clustalw && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    # Install supercronic
    curl -fsSLO "$SUPERCRONIC_URL" && \
    echo "${SUPERCRONIC_SHA1SUM}  ${SUPERCRONIC}" | sha1sum -c - && \
    chmod +x "$SUPERCRONIC" && \
    mv "$SUPERCRONIC" "/usr/local/bin/supercronic"

# Python dependencies (copied and installed first as they change less frequently)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Create cron directories and set up logging (combined directory creation and cron setup)
RUN mkdir -p /var/run/crond /var/log/cron $HOME/cron $HOME/.cache && \
    touch /var/log/cron/cron.log && \
    # Create crontab file for supercronic
    echo "0 * * * * cd $HOME && python -m src.utils.telemetry update_ip_locations >> $HOME/cron/cron.log 2>&1" > $HOME/cron/crontab && \
    echo "*/15 * * * * cd $HOME && curl -X POST http://localhost:8000/api/refresh-telemetry >> $HOME/cron/cron.log 2>&1" >> $HOME/cron/crontab && \
    # Set permissions
    chown -R $USER:$USER $HOME /var/run/crond /var/log/cron

# Copy application code (changes most frequently, so do this last)
COPY ./ ./
RUN chmod +x start-script.sh

# Switch to user
USER $USER

# Expose the application port
EXPOSE 8000

# Start the container
ENTRYPOINT ["./start-script.sh"]