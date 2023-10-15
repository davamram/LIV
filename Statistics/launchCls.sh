#!/bin/bash
for ((value=2000; value<=2400; value+=25))
do
	root -l -q "cls_multiBin.cpp($value)" >> results/cls_etr_$value.txt &
done
