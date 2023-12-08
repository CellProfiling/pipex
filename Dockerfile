FROM python:3.11.6-slim

WORKDIR /opt/pipex

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN apt-get update && apt-get install -y python3-opencv
RUN apt-get install -y libvips

# Temporary install gcc, needed for some of the requirements
RUN apt-get install -y gcc g++

RUN pip install --upgrade pip
COPY ./requirements.txt .
RUN pip install -r requirements.txt

# Remove gcc
RUN apt-get purge -y --auto-remove gcc g++

COPY *.py ./
COPY *.html ./

