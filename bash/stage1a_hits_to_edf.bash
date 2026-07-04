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

python ../py/hits_to_edf.py \
	--hits $hits \
	--$score_or_evalue \
	--truth $truth \
	--lookup ../info/$db.lookup \
	--derived-info ../derived_info/$db.json \
	--algo $algo \
	--reference $db \
	--output $algo.$db.$truth
}

foldseek_hits_cath40=$c/data/big_hits/foldseek.cath40	# q,t,E

foldseek_hits_scop40=$c/data/big_hits/foldseek.scop40	# q,t,E
dali_hits_scop40=$c/data/big_hits/dali.scop40			# q,t,Z
tm_hits_scop40=$c/data/big_hits/tm.scop40				# q,t,TM

reseek28_hits_scop40=$c/data/big_hits/reseek28_sensitive.scop40	# q,t,P

reseek3_kappa_hits_fam_scop40=../hits/reseek_kappa_scop40_fam.hits # q,t,P
reseek3_kappa_hits_fam_cath40=../hits/reseek_kappa_cath40_fam.hits # q,t,P

reseek3_kappa_hits_sf_scop40=../hits/reseek_kappa_scop40_sf.hits # q,t,P
reseek3_kappa_hits_sf_cath40=../hits/reseek_kappa_cath40_sf.hits # q,t,P

reseek3_kappa_hits_fold_scop40=../hits/reseek_kappa_scop40_fold.hits # q,t,P
reseek3_kappa_hits_fold_cath40=../hits/reseek_kappa_cath40_fold.hits # q,t,P

for truth in family superfamily superfamilyx fold
do
	#       algo hits                  fields     s/e     db  truth
	#===============================================================
	run     dali $dali_hits_scop40      1,2,3   score scop40  $truth
	run     dali $dali_hits_scop40      1,2,3   score scop40x $truth

	run     tm   $tm_hits_scop40        1,2,3   score scop40  $truth
	run     tm   $tm_hits_scop40        1,2,3   score scop40x $truth

	run foldseek $foldseek_hits_scop40  1,2,3  evalue scop40  $truth
	run foldseek $foldseek_hits_scop40  1,2,3  evalue scop40x $truth
	run foldseek $foldseek_hits_cath40  1,2,3  evalue cath40  $truth

	run reseek28 $reseek28_hits_scop40  1,2,3  evalue scop40  $truth

done

run reseek3 $reseek3_kappa_hits_fam_scop40  1,2,3  evalue scop40 family
run reseek3 $reseek3_kappa_hits_sf_scop40   1,2,3  evalue scop40 superfamily
run reseek3 $reseek3_kappa_hits_sf_scop40   1,2,3  evalue scop40 superfamilyx
run reseek3 $reseek3_kappa_hits_fold_scop40 1,2,3  evalue scop40 fold

run reseek3 $reseek3_kappa_hits_fam_scop40  1,2,3  evalue scop40x family
run reseek3 $reseek3_kappa_hits_sf_scop40   1,2,3  evalue scop40x superfamily
run reseek3 $reseek3_kappa_hits_sf_scop40   1,2,3  evalue scop40x superfamilyx
run reseek3 $reseek3_kappa_hits_fold_scop40 1,2,3  evalue scop40x fold

run reseek3 $reseek3_kappa_hits_fam_cath40  1,2,3  evalue cath40 family
run reseek3 $reseek3_kappa_hits_sf_cath40   1,2,3  evalue cath40 superfamily
run reseek3 $reseek3_kappa_hits_sf_cath40   1,2,3  evalue cath40 superfamilyx
run reseek3 $reseek3_kappa_hits_fold_cath40 1,2,3  evalue cath40 fold
