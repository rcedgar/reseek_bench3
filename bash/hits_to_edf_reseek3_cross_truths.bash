#!/bin/bash -e

for truth_algo in family superfamily fold
do
	for truth_eval in family superfamily fold
	do
		./hits_to_edf_reseek3_cross_truth.bash $truth_algo $truth_eval
	done
done
