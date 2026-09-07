#!/bin/bash -e

# Convert SCOP/CATH reseek hits to consensus superfamily domain annotations.

ScriptDir=$(cd "$(dirname "$0")" && pwd)
PalmDir=$(dirname "$ScriptDir")
HitDir=${1:-$PalmDir/search_scop_and_cath}
Py=$PalmDir/py/hits_to_domain_annots.py

for prefix in rdrp cdn
do
	for db in scop cath
	do
		hits=$HitDir/${prefix}_${db}.hits
		out=$HitDir/${prefix}_${db}_annots.tsv
		python3 "$Py" --input "$hits" --output "$out"
		echo "wrote $out"
	done
done
