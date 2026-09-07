#!/bin/bash -e

mkdir -p ../annots

cd ../data

python ../py/map_motifs_to_cifs.py \
	rdrp \
	--output ../annots/rdrp_motif_annots.tsv \
	--report ../annots/rdrp_motifs_annots_report.txt

python ../py/map_motifs_to_cifs.py cdn  \
	--output ../annots/cdn_motif_annots.tsv\
	--report ../annots/cdn_motifs_annots_report.txt

python ../py/cifs_to_fasta.py \
	rdrp_cif_chain.files \
	> rdrp_cif.fa

python ../py/cifs_to_fasta.py \
	cdn_cif_chain.files \
	> cdn_cif.fa
