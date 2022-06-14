#!/bin/sh

if [ -f .env ]
then
  export $(cat .env | grep -v '#' | sed 's/#.*//g' | xargs)
fi

bin/python -u pipexGUI.py 2>&1 | tee log.txt

