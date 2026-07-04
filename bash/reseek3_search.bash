#!/bin/bash -e

mkdir -p ../hits
mkdir -p ../time

for ref in scop40 cath40
do
	db=../reseek_dbs/$ref.bcb
	for mode in kappa # all
	do
		for truth in fam sf fold
		do
			name=reseek_${mode}_${ref}_$truth
			hits=../hits/$name.hits

			/bin/time -vo ../time/$name \
			reseek \
				-flat_search_$mode $db \
				-stats $truth \
				-db $db \
				-output $hits \
				-log ../hits/$name.search.log
		done
	done
done
