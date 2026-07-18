#!/bin/bash -e

big_outdir=../big_invidx_test
outdir=../invidx_test

mkdir -p $outdir $big_outdir

bcb=../big_reseek_dbs/scop40x.bcb
kdx=$big_outdir/scop40x.kdx
opts="-columns query+target+pvalue -stats sf -fast"

cd $outdir

reseek -createindex $bcb -output $kdx -log createindex.log

reseek -idx_search_kappa $bcb  -db $kdx -input2 $bcb $opts -log kdx.log -output kdx.hits 
reseek -flat_search_kappa $bcb -db $bcb              $opts -log qdx.log -output qdx.hits

reseek -fast_bench_hits kdx.hits -truth sf -log bench_kdx.log
reseek -fast_bench_hits qdx.hits -truth sf -log bench_qdx.log

grep Sum3= bench_?dx.log

grep -i "ax mem" ?dx.log
grep -i "lapsed" ?dx.log
