#!/bin/bash -e

ref=scop40

db=../big_foldseek_dbs/$ref/$ref
tmpdir=../tmp_lol
mkdir -p ../tmp_lol


foldseek easy-search \
    $db \
    $db \
    ../big_hits/lol.scop40.hits \
    $tmpdir \
    --alignment-type 3 \
    --format-mode 4 \
    --max-seqs 2000 \
    --format-output "query,target,lddt,qtmscore,alntmscore,evalue,bits,fident,qstart,qend,tstart,tend,cigar" \
    --threads 64
