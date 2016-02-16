#!/bin/sh
mkdir pico
i="0"
while [ $i -lt 10 ]
do
	tar -xvf files/Files_$1_$i.tar
	cp files/Files_$1_$i/picodst/* pico
	rm -rf files/Files_$1_$i
	i=$[$i+1]
done

