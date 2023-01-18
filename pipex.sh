#!/bin/sh

if [ -f .env ]
then
  export $(cat .env | grep -v '#' | sed 's/#.*//g' | xargs)
fi

if test -f "./bin/python"
then
  ./bin/python -u pipexGUI.py 2>&1 | tee log.txt
else
  python -u pipexGUI.py 2>&1 | tee log.txt
fi
