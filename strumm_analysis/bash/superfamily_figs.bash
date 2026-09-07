#!/bin/bash -e

# Generate SCOP/CATH superfamily domain figures for rdrp and cdn annots.
# Intended to run under WSL with conda env null_model.

ScriptDir=$(cd "$(dirname "$0")" && pwd)
PalmDir=$(dirname "$ScriptDir")
AnnotDir=$PalmDir/search_scop_and_cath
OutDir=$PalmDir/figs
Py=$PalmDir/py/superfamily_fig.py

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
	# shellcheck source=/dev/null
	source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
	# shellcheck source=/dev/null
	source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
	eval "$(conda shell.bash hook)"
fi
conda activate null_model

Reseek=$PalmDir/binaries/reseek
if [ ! -x "$Reseek" ]; then
	Reseek=$(command -v reseek || true)
fi
if [ -z "$Reseek" ]; then
	echo "*ERROR* reseek not found" >&2
	exit 1
fi

PalmcoreFa=$PalmDir/data/rdrp_abc_palmcore.fa
"$Reseek" -convert "$PalmDir/data/rdrp_abc_palmcore.bcb" -fasta "$PalmcoreFa"

mkdir -p "$OutDir"

# Palmcore first so annotated rdrp tips keep search coordinates; full-length
# fills in order entries (e.g. decoys) absent from palmcore.
python3 "$Py" \
	--domains "$AnnotDir/rdrp_cath_annots.tsv" \
	--motifs "$PalmDir/data/rdrp_motifs.tsv" \
	--fasta "$PalmcoreFa" \
	--fasta "$PalmDir/data/rdrp.fa" \
	--motif-style abc \
	--order "$PalmDir/data/rdrp_figure_order.txt" \
	--sf-names "$PalmDir/data/cath-b-newest-names" \
	--cath \
	--view-lo 1 \
	--output "$OutDir/rdrp_cath_sf.pdf"

python3 "$Py" \
	--domains "$AnnotDir/rdrp_scop_annots.tsv" \
	--motifs "$PalmDir/data/rdrp_motifs.tsv" \
	--fasta "$PalmcoreFa" \
	--fasta "$PalmDir/data/rdrp.fa" \
	--motif-style abc \
	--order "$PalmDir/data/rdrp_figure_order.txt" \
	--sf-names "$PalmDir/data/scop_superfamily_names.txt" \
	--view-lo 1 \
	--output "$OutDir/rdrp_scop_sf.pdf"

python3 "$Py" \
	--domains "$AnnotDir/cdn_cath_annots.tsv" \
	--motifs "$PalmDir/data/cdn_motifs.tsv" \
	--fasta "$PalmDir/data/cdn.fa" \
	--motif-style iii \
	--order "$PalmDir/data/cdn_figure_order.txt" \
	--sf-names "$PalmDir/data/cath-b-newest-names" \
	--cath \
	--view-lo 1 \
	--output "$OutDir/cdn_cath_sf.pdf"

python3 "$Py" \
	--domains "$AnnotDir/cdn_scop_annots.tsv" \
	--motifs "$PalmDir/data/cdn_motifs.tsv" \
	--fasta "$PalmDir/data/cdn.fa" \
	--motif-style iii \
	--order "$PalmDir/data/cdn_figure_order.txt" \
	--sf-names "$PalmDir/data/scop_superfamily_names.txt" \
	--view-lo 1 \
	--output "$OutDir/cdn_scop_sf.pdf"

echo "done: $OutDir/{rdrp,cdn}_{cath,scop}_sf.pdf"
