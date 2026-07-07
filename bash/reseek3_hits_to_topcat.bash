#!/bin/bash -e

# Fast stage1b: reuse binary cache, C++ hits_to_topcat --bin.
# Original stage1b_hits_to_topcat.bash (Python) is unchanged.

CPP_BIN=../cpp/bin
HITS_TO_BIN=$CPP_BIN/hits_to_bin
HITS_TO_TOPCAT=$CPP_BIN/hits_to_topcat
CACHE_DIR=../cache/bin

if [ ! -x "$HITS_TO_BIN" ] || [ ! -x "$HITS_TO_TOPCAT" ]; then
	echo "Build C++ tools first: cd cpp && make"
	exit 1
fi

bin_cache_path () {
	hits=$1
	db=$2
	base=$(basename "$hits")
	base=${base%%.*}
	echo "$CACHE_DIR/${base}.${db}.bin"
}

ensure_bin () {
	hits=$1
	db=$2
	fields=$3
	binpath=$(bin_cache_path "$hits" "$db")
	mkdir -p "$CACHE_DIR"
	if [ ! -f "$binpath" ] || [ "$hits" -nt "$binpath" ]; then
		echo "hits_to_bin: $hits -> $binpath" >&2
		"$HITS_TO_BIN" \
			--hits "$hits" \
			--lookup "../info/$db.lookup" \
			--fields "$fields" \
			--output "$binpath"
	fi
	echo "$binpath"
}

run_one_truth () {
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

binpath=$(ensure_bin "$hits" "$db" "$fields")

"$HITS_TO_TOPCAT" \
	--bin "$binpath" \
	--$score_or_evalue \
	--truth "$truth" \
	--lookup "../info/$db.lookup" \
	--derived-info "../derived_info/$db.json" \
	--algo "$algo" \
	--reference "$db" \
	--output "$algo.$db.$truth"
}

run_all_truths () {
algo=$1
hits=$2
fields=$3
score_or_evalue=$4
db=$5

mkdir -p ../edf
cd ../edf

binpath=$(ensure_bin "$hits" "$db" "$fields")

"$HITS_TO_TOPCAT" \
	--bin "$binpath" \
	--$score_or_evalue \
	--all-truths \
	--lookup "../info/$db.lookup" \
	--derived-info "../derived_info/$db.json" \
	--algo "$algo" \
	--reference "$db" \
	--output "$algo.$db"
}

foldseek_hits_cath40=$c/data/big_hits/foldseek.cath40.hits	# q,t,E

foldseek_hits_scop40=$c/data/big_hits/foldseek.scop40.hits	# q,t,E
dali_hits_scop40=$c/data/big_hits/dali.scop40			# q,t,Z
tm_hits_scop40=$c/data/big_hits/tm.scop40				# q,t,TM

reseek28_hits_scop40=$c/data/big_hits/reseek28_sensitive.scop40	# q,t,P

reseek3_fast_hits_sf_scop40=../big_hits/reseek_fast_scop40_sf.hits # q,t,P
reseek3_fast_hits_sf_cath40=../big_hits/reseek_fast_cath40_sf.hits # q,t,P

reseek3_fast_hits_fold_scop40=../big_hits/reseek_fast_scop40_fold.hits # q,t,P
reseek3_fast_hits_fold_cath40=../big_hits/reseek_fast_cath40_fold.hits # q,t,P

reseek3_sensitive_hits_sf_scop40=../big_hits/reseek_sensitive_scop40_sf.hits # q,t,P
reseek3_sensitive_hits_sf_cath40=../big_hits/reseek_sensitive_cath40_sf.hits # q,t,P

reseek3_sensitive_hits_fold_scop40=../big_hits/reseek_sensitive_scop40_fold.hits # q,t,P
reseek3_sensitive_hits_fold_cath40=../big_hits/reseek_sensitive_cath40_fold.hits # q,t,P



#       algo hits                  fields     s/e     db
#===============================================================
run_one_truth reseek3 $reseek3_fast_hits_sf_scop40   1,2,3  evalue scop40 topsf
run_one_truth reseek3 $reseek3_fast_hits_fold_scop40 1,2,3  evalue scop40 topfold

run_one_truth reseek3 $reseek3_fast_hits_sf_scop40   1,2,3  evalue scop40x topsf
run_one_truth reseek3 $reseek3_fast_hits_fold_scop40 1,2,3  evalue scop40x topfold

run_one_truth reseek3 $reseek3_fast_hits_sf_cath40 1,2,3  evalue cath40 topfold
run_one_truth reseek3 $reseek3_fast_hits_fold_cath40 1,2,3  evalue cath40 topfold

run_one_truth reseek3 $reseek3_sensitive_hits_sf_scop40   1,2,3  evalue scop40 topsf
run_one_truth reseek3 $reseek3_sensitive_hits_fold_scop40 1,2,3  evalue scop40 topfold

run_one_truth reseek3 $reseek3_sensitive_hits_sf_scop40   1,2,3  evalue scop40x topsf
run_one_truth reseek3 $reseek3_sensitive_hits_fold_scop40 1,2,3  evalue scop40x topfold

run_one_truth reseek3 $reseek3_sensitive_hits_sf_cath40 1,2,3  evalue cath40 topfold
run_one_truth reseek3 $reseek3_sensitive_hits_fold_cath40 1,2,3  evalue cath40 topfold



