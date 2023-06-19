#!/bin/bash

clear
FILENAME="PATCON"
for i in {1..2..1}
do
    date
    echo "STARTING THE SIMULATION NUMBER $i"
    sleep 1
    echo "COPYING FILE"
    cp ./input/$FILENAME$i.txt ./input/$FILENAME.txt
    sleep 0.5
    echo "RUNNING THE SIMULATION"
    echo "START!"
    python src/main.py
    echo "SIMULATION NUMBER $i FINISHED"
done
echo "SIMULATION DONE"