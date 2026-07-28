
c:/src/reseek/github_releases/reseek-v3.01-win64.exe \
	-flat_search_kappa ../big_reseek_dbs/cath40.bcb \
	-sensitive \
	-stats fold \
	-db ../big_reseek_dbs/cath40.bcb \
	-output ../big_hits/reseek3_sensitive_fold.cath40.hits \
	-columns query+target+pvalue \
	-log ../logs/reseek_v3.01_win64_fast_cath40_fold.log
