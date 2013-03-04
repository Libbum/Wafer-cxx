#!/bin/sh

./cleandatafiles.sh
cp params.txt data
date
nohup nice make run
date
