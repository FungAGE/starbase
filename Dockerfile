# Select base image
FROM python:3.9
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"

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
    apt-get install -y wget \
    build-essential \
    cmake \
    zlib1g-dev \
    liblmdb-dev \
    libncurses5-dev \
    libboost-all-dev \
    ncbi-blast+ \
    hmmer \
    clustalw && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Python dependencies separately (cache this layer)
RUN pip install --no-cache-dir -r requirements.txt

# # Install BLAST
# RUN wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-src.tar.gz && \
#     tar -zxzf ncbi-blast-2.9.0+-src.tar.gz && \
#     cd ncbi-blast-2.9.0+-src/c++ && \
#     ./configure --prefix=/usr/local/BLAST2.9 --without-debug --without-exe --without-boost --without-gui && \
#     make -j8 && \
#     make install && \
#     cd ../../.. && \
#     rm -rf ncbi-blast-2.9.0+-src ncbi-blast-2.9.0+-src.tar.gz

# Install DIAMOND
RUN wget -q https://github.com/bbuchfink/diamond/archive/v2.1.9.tar.gz && \
    tar xzf v2.1.9.tar.gz && \
    mkdir diamond-2.1.9/bin && \
    cd diamond-2.1.9/bin && \
    cmake -DBLAST_INCLUDE_DIR=/usr/local/BLAST2.9/include/ncbi-tools++ .. && \
    make -j8 && \
    make install && \
    cd ../../.. && \
    rm -rf diamond-2.1.9 diamond-2.1.9.tar.gz

# Set the PATH environment variable to include the BLAST and DIAMOND binaries
ENV PATH="/usr/local/BLAST2.9/bin:${PATH}"

# Copy the rest of the code
COPY ./ ./

# unzip db
RUN gunzip -c src/data/db.tar.gz | tar -xf - -C src/data

# build blast dbs from sql table
RUN python /src/utils/blastdb.py

# Change permissions
RUN chmod +x start-script.sh && chown -R $USER:$USER $HOME

USER $USER

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]
