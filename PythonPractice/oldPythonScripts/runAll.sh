#!/bin/sh

echo "Hello World!"
echo "knfdknfg"

yearFrom=2015
yearTo=2022

name1='forcing.'
name2='.log'
for (( i=$yearFrom; i<=$yearTo; i++ ))
do
echo $i
FILE="$name1$i$name2"
if [ ! -f $FILE ]
then
    sbatch -J met --export=ALL,yearToRun=${i} submitGetData.bsub
fi
done
