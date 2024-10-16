# Select base image
FROM python:3.9
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="starbase is a database and toolkit for exploring of large transposable elements in  Fungi"

# Create user name and home directory variables
ENV USER=starbase
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory
WORKDIR $HOME/

# Copy only the requirements.txt first for dependency installation
COPY requirements.txt .

# Update system and install system dependencies first
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y ncbi-blast+ hmmer clustalw && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python dependencies separately (cache this layer)
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the code
COPY ./ ./

# unzip db
RUN cd src/data && gunzip -c db.tar.gz | tar -xf - && ls -l db/ && chmod 644 db/starbase.sqlite && cd ../../
# to create the db.tar.gz 
# cd src/data/ && tar -cf db.tar db/* && gzip -f db.tar

# build blast dbs from sql table
RUN python src/utils/blastdb.py

# Change permissions
RUN chmod +x start-script.sh && chown -R $USER:$USER $HOME

USER $USER

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]
