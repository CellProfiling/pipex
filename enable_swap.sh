#!/bin/sh

SWAP_REQ=`expr $1`
SWAP_REQ=`expr $SWAP_REQ \* 1024 \* 1024 \* 1024`

SWAP_SIZE=$(free | grep "Swap" | awk '{print $2}')
if [ $SWAP_REQ -gt $SWAP_SIZE ] 
then
    dd if=/dev/zero of=./swapfile count=$1K bs=1M
    mkswap ./swapfile
    sudo chown root:root ./swapfile
    sudo chmod 600 ./swapfile
    sudo swapon ./swapfile
fi

