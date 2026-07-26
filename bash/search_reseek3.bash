#!/bin/bash -e

#################################
rm -rf ../big_hits/reseek3_*.hits
rm -rf ../logs/reseek3_*.search.log
###################################

mkdir -p ../big_hits
mkdir -p ../time

for ref in scop40 cath40
do
	db=../big_reseek_dbs/$ref.bcb
	for mode in fast sensitive verysensitive
	do
		for truth in family superfamily fold
		do
			name=reseek3_${mode}_${truth}.${ref}
			hits=../big_hits/$name.hits

			cmd="reseek \
					-flat_search_kappa $db \
					-$mode \
					-stats $truth \
					-db $db \
					-output $hits \
					-columns query+target+pvalue \
					-log ../logs/$name.search.log"

			if [ -x /bin/time ] ; then
				/bin/time -v -o ../time/$name.search.time $cmd
			else
				$cmd
			fi
		done
	done
done
echo ===========================
echo SECONDS=$SECONDS \
	| tee ../logs/$name.search.seconds
echo MINUTES=$(($SECONDS/60))
echo ===========================
