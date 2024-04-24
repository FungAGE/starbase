# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.9
MAINTAINER Adrian Forsythe <adrian.e.forsythe@gmail.com>

# Create user name and home directory variables. 
# The variables are later used as $USER and $HOME. 
ENV USER=starbase
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/

# Copy code and start script (this will place the files in home/username/)
COPY ./ ./

# Update system and install dependencies.
RUN apt-get update && apt-get upgrade -y && apt-get install ncbi-blast+ hmmer -y && apt-get clean && rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir -r $HOME/requirements.txt && \
    chmod +x start-script.sh && \
    # mkdir database_folder/ && \
    chown -R $USER:$USER $HOME

USER $USER

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]
