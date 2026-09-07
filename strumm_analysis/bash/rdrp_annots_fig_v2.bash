#!/bin/bash -e

# Intended to run under WSL (matplotlib via conda env null_model for user bob)

ScriptDir=$(cd "$(dirname "$0")" && pwd)
PalmDir=$(dirname "$ScriptDir")
OutDir=$PalmDir/figs
PyDir=$PalmDir/py

mkdir -p "$OutDir"

# Prefer bob's conda env (local .venv* removed).
NullPy=/home/bob/miniconda3/envs/null_model/bin/python
if [ -x "$NullPy" ]; then
	Py=$NullPy
elif command -v conda >/dev/null 2>&1; then
	# shellcheck source=/dev/null
	source "$(conda info --base)/etc/profile.d/conda.sh"
	conda activate null_model
	Py=python
else
	Py=python3
fi

"$Py" "$PyDir/make_annots_fig_v2.py" --family rdrp \
	--dbtypes SCOP,Pfam \
	--output "$OutDir/rdrp_annots_v2.pdf,$OutDir/rdrp_annots_v2.png"

ls -lh "$OutDir/rdrp_annots_v2.pdf" "$OutDir/rdrp_annots_v2.png"
