#!/bin/bash

# p restart iter len
for file in bacp8-2.txt toy.txt; do
	for i in 1 2 0; do
		for j in 10 50; do
			for k in 200 1000; do
				for l in 1 2 3; do
					while read line; do
  				./bacp $file $i $j $k $l $line >> results/results_${file}_${i}_${j}_${k}_${l}.txt
					done <seeds.txt
				done
			done
		done
	done
done

