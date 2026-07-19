#!/bin/bash -e

big_outdir=../big_invidx_test
outdir=../invidx_test

##############
rm -rf $outdir
rm -rf $big_outdir
##################

mkdir -p $outdir $big_outdir

bcb=../big_reseek_dbs/scop40x.bcb
kdx=$big_outdir/scop40x.kdx
opts="-idx_search_kappa $bcb -db $bcb -columns query+target+pvalue -stats sf -fast -threads 64"

cd $outdir

reseek -createindex $bcb -output $kdx -log createindex.log

reseek $opts -log search_max.log   -output shards.hits   -unique_kmer -max_seqs 1500 -db_shards 64
reseek $opts -log search_max.log   -output max.hits   -unique_kmer -max_seqs 1500
reseek $opts -log search_ukmer.log -output ukmer.hits -unique_kmer
reseek $opts -log search_kdx.log   -output kdx.hits   -kdx $kdx
reseek $opts -log search_nokdx.log -output nokdx.hits
reseek $opts -log search_qdx.log   -output qdx.hits

reseek -fast_bench_hits shards.hits -truth sf -log bench_shards.log
reseek -fast_bench_hits max.hits    -truth sf -log bench_max.log
reseek -fast_bench_hits ukmer.hits  -truth sf -log bench_ukmer.log
reseek -fast_bench_hits kdx.hits    -truth sf -log bench_kdx.log
reseek -fast_bench_hits nokdx.hits  -truth sf -log bench_nokdx.log
reseek -fast_bench_hits qdx.hits    -truth sf -log bench_qdx.log

grep Sum3= bench_*.log

grep -i "ax mem" search_*.log
grep -i "lapsed" search_*.log
