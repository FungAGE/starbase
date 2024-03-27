FROM python:3.9

RUN apt-get update && apt-get upgrade -y && apt-get clean && rm -rf /var/lib/apt/lists/*
  
COPY ./requirements.txt /requirements.txt

RUN pip3 install -r /requirements.txt

RUN mkdir /starbase
WORKDIR /starbase
COPY ./ ./

ENV GUNICORN_CMD_ARGS "--bind=0.0.0.0:8000 --workers=2 --thread=4 --worker-class=gthread --forwarded-allow-ips='*' --access-logfile -"

CMD [ "gunicorn", "src.app:server"]