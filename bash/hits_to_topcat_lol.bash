#!/bin/bash -e

#################################################################################
# https://github.com/steineggerlab/foldseek#output-search [39b02b9f]
#
# The e-value field #11 reports the normalized LoLscore (<= 1)
# score/bits field #12 reports the estimated homology probability scaled by 10000
#	(i.e. probability × 10000, range 0–10000).
#################################################################################

mkdir -p ../topcat

for mode in prefilter exhaustive
do
	for ref in scop40 scop40x
	do
		for truth in topfold topsf
		do
			python ../py/hits_to_topcat.py \
				--hits ../big_hits/lol_${mode}.scop40.allfields.hits \
				--algo lol_score_prefilter \
				--lookup ../data/$ref.lookup \
				--derived-info ../derived_info/$ref.json \
				--truth $truth \
				--fields 1,2,11 \
				--score \
				--output ../topcat/lol_score_${mode}.$ref.$truth

			python ../py/hits_to_topcat.py \
				--hits ../big_hits/lol_${mode}.scop40.allfields.hits \
				--algo lol_prob_prefilter \
				--lookup ../data/$ref.lookup \
				--derived-info ../derived_info/$ref.json \
				--truth $truth \
				--fields 1,2,12 \
				--score \
				--output ../topcat/lol_prob_${mode}.$ref.$truth
		done
	done
done
