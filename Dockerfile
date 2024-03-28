FROM python:3.9

RUN apt-get update && apt-get upgrade -y && apt-get clean && rm -rf /var/lib/apt/lists/*
COPY ./requirements.txt /requirements.txt

RUN pip3 install -r /requirements.txt

RUN mkdir /home/starbase/
RUN mkdir /home/project-vol/
WORKDIR /home/starbase/
COPY ./ ./

COPY ./start-script.sh ./src/
RUN chmod +x ./src/start-script.sh

EXPOSE 8000

ENTRYPOINT ["./src/start-script.sh"]