# Select base image
FROM python:3.9
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="starbase is a database and toolkit for exploring large transposable elements in Fungi"

# Create variables for user name, home directory, and placeholders 
ENV USER=starbase
ENV HOME=/home/$USER
ENV DB_USER=""
ENV DB_PASSWORD=""
ENV DB_HOST=""
ENV DB_PORT="3306"
ENV DB_NAME=""
ENV TS_AUTH_KEY=""

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory
WORKDIR $HOME/

# Copy only the requirements.txt first for dependency installation
COPY requirements.txt .

# Update system and install system dependencies first
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y curl iptables ncbi-blast+ hmmer clustalw && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python dependencies separately (cache this layer)
RUN pip install --no-cache-dir -r requirements.txt

# Install Tailscale
RUN curl -fsSL https://tailscale.com/install.sh | sh

# Copy the rest of the code
COPY ./ ./

# Run the Tailscale login script
RUN chmod +x tailscale.sh && chmod +x start-script.sh && sh tailscale.sh

# Change permissions for user
RUN chown -R $USER:$USER $HOME

# Switch to user
USER $USER

# Expose the application port
EXPOSE 8000

# Start the container by initializing Tailscale and running the main app script
ENTRYPOINT ["./start-script.sh"]
