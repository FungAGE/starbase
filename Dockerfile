# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.9

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

RUN mkdir /home/project-vol/

# Update system and install dependencies.
RUN apt-get update && apt-get upgrade -y && apt-get clean && rm -rf /var/lib/apt/lists/* && \
    chmod +x start-script.sh && \
    chown -R $USER:$USER $HOME

USER $USER

RUN pip install --no-cache-dir -r $HOME/requirements.txt \

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]