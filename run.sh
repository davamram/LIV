#!/bin/bash
for ((value=2000; value<=2400; value+=25))
do
	./launch.sh $value
	./plot.sh $value
done
