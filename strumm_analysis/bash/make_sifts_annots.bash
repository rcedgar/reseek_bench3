#!/bin/bash -e

mkdir -p ../annots
cd ../data

for fam in cdn rdrp
do
	python ../py/cif_make_sifts_annotations.py \
		--input ${fam}_cif_chain.files \
		--output ../annots/${fam}_sifts_annots.tsv \
		--report ../annots/${fam}_sifts_annots_report.txt
done
