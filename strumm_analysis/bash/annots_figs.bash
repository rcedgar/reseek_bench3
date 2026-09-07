#!/bin/bash -e

# Reformat SIFTS annots and draw rdrp/cdn annotation figures.
# Intended to run under WSL (matplotlib via .venv_wsl or system python3).

ScriptDir=$(cd "$(dirname "$0")" && pwd)
PalmDir=$(dirname "$ScriptDir")
OutDir=$PalmDir/figs
PyDir=$PalmDir/py

mkdir -p "$OutDir"

Py=python3
if [ -x "$PalmDir/.venv_wsl/bin/python" ]; then
	Py=$PalmDir/.venv_wsl/bin/python
fi

echo "reformat SIFTS fig annots..."
"$Py" "$PyDir/reformat_sifts_annots.py" --family rdrp
"$Py" "$PyDir/reformat_sifts_annots.py" --family cdn

echo "draw figures..."
"$Py" "$PyDir/make_annots_fig.py" --family rdrp \
	--dbtypes SCOP,Pfam \
	--output "$OutDir/rdrp_annots.pdf,$OutDir/rdrp_annots.png"
"$Py" "$PyDir/make_annots_fig.py" --family cdn \
	--output "$OutDir/cdn_annots.pdf,$OutDir/cdn_annots.png"

echo "done: $OutDir/{rdrp,cdn}_annots.pdf/.png"
