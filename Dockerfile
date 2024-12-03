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
    apt-get install -y curl iptables ncbi-blast+ hmmer clustalw && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY ./ ./

# Make script executable
RUN chmod +x start-script.sh

# Run precomputation after code is available
RUN python3 -c "from src.components.sql_manager import precompute_all; precompute_all()"

# Change permissions for user
RUN chown -R $USER:$USER $HOME

# Switch to user
USER $USER

# Expose the application port
EXPOSE 8000

# Start the container
ENTRYPOINT ["./start-script.sh"]