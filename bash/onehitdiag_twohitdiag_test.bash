#!/bin/bash -e
 
# reseek [d2e70a8]
 
# SEPQ0.1=0.298 SEPQ1=0.407 SEPQ10=0.462 Sum3=1.667 sf Kappa filter 13 secs	9.80Gb	01:16	default	
# SEPQ0.1=0.298 SEPQ1=0.407 SEPQ10=0.462 Sum3=1.667 sf Kappa filter 19 secs	9.68Gb	01:19	onehitdiag
# SEPQ0.1=0.287 SEPQ1=0.383 SEPQ10=0.424 Sum3=1.573 sf Kappa filter  7 secs	9.70Gb	01:02	twohitdiag
 
# PDB vs. 100
# Kappa filter 39 secs	4.91Gb	2:09	default
# Kappa filter 71 secs	4.12Gb	3:19	onehitdiag
# Kappa filter 22 secs	4.52Gb	1:48	twohitdiag
 
# AFDB50 vs. 100
# Kappa filter 252 secs	10.2Gb	4:48	default
# Kappa filter 472 secs	11.3Gb	08:14	onehitdiag
# Kappa filter 121 secs	10.2Gb	02:46	twohitdiag

#######################################################
# SCOP40x
#######################################################
q=../data/scop40x.bcb
db=/data/scop40x.bcb

reseek \
	-flat_search_kappa $q \
	-db $db \
	-fast \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-log default.log \
	-output default.hits

reseek \
	-flat_search_kappa $q \
	-db $db \
	-fast \
	-onehitdiag
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-log onehitdiag.log \
	-output onehitdiag.hits

reseek \
	-flat_search_kappa $q \
	-db $db \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-fast \
	-twohitdiag \
	-log twohitdiag.log \
	-output twohitdiag.hits

#######################################################
# Bench SCOP40x
#######################################################
db=../data/scop40x.bcb
lookup=../data/scop40x.lookup
truth=sf

for name in default onehitdiag twohitdiag
do
	reseek \
		-fast_bench_hits $name.hits \
		-lookup $lookup \
		-truth $truth \
		-log $name.bench.log
done

#######################################################
# PDB
#######################################################
q=c:/src/reseek_bench3/big_reseek_dbs/pdb_subset1000.bcb
db=c:/data/pdb/reseek_db/pdb.bcb

reseek \
	-flat_search_kappa $q \
	-db $db \
	-fast \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-log pdb_default.log \
	-output pdb_default.hits

reseek \
	-flat_search_kappa $q \
	-db $db \
	-fast \
	-onehitdiag \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-log pdb_onehitdiag.log \
	-output pdb_onehitdiag.hits

reseek \
	-flat_search_kappa $q \
	-db $db \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-fast \
	-twohitdiag \
	-log pdb_twohitdiag.log \
	-output pdb_twohitdiag.hits

#######################################################
# AFDB50
#######################################################
q=c:/src/reseek_bench3/big_reseek_dbs/pdb_subset100.bcb
db=c:/data/afdb/afdb50.bcb

reseek \
	-flat_search_kappa $q \
	-db $db \
	-fast \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-log afdb50_default.log \
	-output afdb50_default.hits

reseek \
	-flat_search_kappa $q \
	-db $db \
	-fast \
	-onehitdiag \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-log afdb50_onehitdiag.log \
	-output afdb50_onehitdiag.hits

reseek \
	-flat_search_kappa $q \
	-db $db \
	-stats sf \
	-threads 90 \
	-columns query+target+pvalue \
	-fast \
	-twohitdiag \
	-log afdb50_twohitdiag.log \
	-output afdb50_twohitdiag.hits
