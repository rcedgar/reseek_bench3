#!/bin/bash -e

run () {
algo=$1
hits=$2
fields=$3
score_or_evalue=$4
db=$5
truth=$6

if [ -z "$truth" ] ; then
	echo Missing arg
	exit 1
fi

mkdir -p ../edf
cd ../edf

python ../py/hits_to_topcat.py \
	--hits $hits \
	--$score_or_evalue \
	--truth $truth \
	--lookup ../info/$db.lookup \
	--derived-info ../derived_info/$db.json \
	--algo $algo \
	--output $algo.$db.$truth
}

foldseek_hits_cath40=$c/data/big_hits/foldseek.cath40.hits	# q,t,E

foldseek_hits_scop40=$c/data/big_hits/foldseek.scop40.hits	# q,t,E
dali_hits_scop40=$c/data/big_hits/dali.scop40			# q,t,Z
tm_hits_scop40=$c/data/big_hits/tm.scop40				# q,t,TM

reseek28_hits_scop40=$c/data/big_hits/reseek28_sensitive.scop40	# q,t,P

reseek3_kappa_hits_fam_scop40=../big_hits/reseek_kappa_scop40_fam.hits # q,t,P
reseek3_kappa_hits_fam_cath40=../big_hits/reseek_kappa_cath40_fam.hits # q,t,P

reseek3_kappa_hits_sf_scop40=../big_hits/reseek_kappa_scop40_sf.hits # q,t,P
reseek3_kappa_hits_sf_cath40=../big_hits/reseek_kappa_cath40_sf.hits # q,t,P

reseek3_kappa_hits_fold_scop40=../big_hits/reseek_kappa_scop40_fold.hits # q,t,P
reseek3_kappa_hits_fold_cath40=../big_hits/reseek_kappa_cath40_fold.hits # q,t,P

for truth in topsf topfold
do
	#       algo hits                  fields     s/e     db  truth
	#===============================================================
	run     dali $dali_hits_scop40      1,2,3   score scop40  $truth
	run     dali $dali_hits_scop40      1,2,3   score scop40x $truth

	run     tm   $tm_hits_scop40        1,2,3   score scop40  $truth
	run     tm   $tm_hits_scop40        1,2,3   score scop40x $truth

	run foldseek $foldseek_hits_scop40  1,2,3  evalue scop40  $truth
	run foldseek $foldseek_hits_scop40  1,2,3  evalue scop40x $truth
	if [ "$truth" = topfold ] ; then
		run foldseek $foldseek_hits_cath40  1,2,3  evalue cath40 topfold
	fi

	run reseek28 $reseek28_hits_scop40  1,2,3  evalue scop40  $truth
done

run reseek3 $reseek3_kappa_hits_sf_scop40   1,2,3  evalue scop40 topsf
run reseek3 $reseek3_kappa_hits_fold_scop40 1,2,3  evalue scop40 topfold

run reseek3 $reseek3_kappa_hits_sf_scop40   1,2,3  evalue scop40x topsf
run reseek3 $reseek3_kappa_hits_fold_scop40 1,2,3  evalue scop40x topfold

run reseek3 $reseek3_kappa_hits_fold_cath40 1,2,3  evalue cath40 topfold
