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

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables ncbi-blast+ hmmer clustalw cron && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY ./ ./

# Make script executable
RUN chmod +x start-script.sh

# Run precomputation after code is available
# RUN python3 -c "from src.components.sql_manager import precompute_all; precompute_all()"

# Update cron job to refresh both IP locations and telemetry data
RUN echo "0 * * * * cd $HOME && python -m src.utils.telemetry update_ip_locations >> /var/log/cron.log 2>&1" > /etc/cron.d/telemetry-cron && \
    echo "*/15 * * * * cd $HOME && curl -X POST http://localhost:8000/api/refresh-telemetry >> /var/log/cron.log 2>&1" >> /etc/cron.d/telemetry-cron

# Give execution rights on the cron job
RUN chmod 0644 /etc/cron.d/telemetry-cron

# Apply the cron job
RUN crontab /etc/cron.d/telemetry-cron

# Create the log file with appropriate permissions
RUN touch /var/log/cron.log && \
    chmod 666 /var/log/cron.log

# Change permissions for user
RUN chown -R $USER:$USER $HOME

# Switch to user
USER $USER

# Expose the application port
EXPOSE 8000

# Start the container
ENTRYPOINT ["./start-script.sh"]