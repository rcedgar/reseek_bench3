#!/bin/bash -e

# Draw rdrp CATH-only annotation figure (overlap-aware sub-levels).
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

echo "draw rdrp CATH figure..."
"$Py" "$PyDir/make_rdrp_cath_fig.py" \
	--output "$OutDir/rdrp_cath_annots.pdf,$OutDir/rdrp_cath_annots.png"

echo "done: $OutDir/rdrp_cath_annots.pdf/.png"
