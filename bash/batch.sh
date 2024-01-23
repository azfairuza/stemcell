#!/bin/bash

clear
FILENAME="CELCON"
FILENAME2="PATCON"
for i in {1..3..1}
do
    date
    echo "STARTING THE SIMULATION NUMBER $i"
    sleep 1
    echo "COPYING FILE"
    cp ./input/$FILENAME$i.txt ./input/$FILENAME.txt
    cp ./input/$FILENAME2$i.txt ./input/$FILENAME2.txt
    sleep 0.5
    echo "RUNNING THE SIMULATION"
    echo "START!"
    python src/main.py
    echo "SIMULATION NUMBER $i FINISHED"
done
echo "SIMULATION DONE"