#!/bin/bash -e

# Combined Pfam + SCOP + CATH domain figures.
# Intended to run under WSL with conda env null_model.

ScriptDir=$(cd "$(dirname "$0")" && pwd)
PalmDir=$(dirname "$ScriptDir")
AnnotDir=$PalmDir/search_scop_and_cath
OutDir=$PalmDir/figs
Py=$PalmDir/py/superfamily_fig2.py
# Pfam annotations (full-length coords; mapped onto tip FASTA when needed).
PfamDir=${PFAM_DIR:-/mnt/c/src/palmfinder/ref_pfam}
CdnPfamDir=${CDN_PFAM_DIR:-/mnt/c/src/mine_gfp_nanoluc/cdntase2/pfam_ref}

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

python3 "$Py" \
	--scop-domains "$AnnotDir/rdrp_scop_annots.tsv" \
	--cath-domains "$AnnotDir/rdrp_cath_annots.tsv" \
	--pfam-domains "$PfamDir/domains.tsv" \
	--pfam-names "$PfamDir/pfams.txt" \
	--motifs "$PalmDir/data/rdrp_motifs.tsv" \
	--fasta "$PalmcoreFa" \
	--fasta "$PalmDir/data/rdrp.fa" \
	--motif-style abc \
	--order "$PalmDir/data/rdrp_figure_order.txt" \
	--scop-names "$PalmDir/data/scop_superfamily_names.txt" \
	--cath-names "$PalmDir/data/cath-b-newest-names" \
	--view-lo 1 \
	--output "$OutDir/rdrp_cath_sf2.pdf"

python3 "$Py" \
	--scop-domains "$AnnotDir/cdn_scop_annots.tsv" \
	--cath-domains "$AnnotDir/cdn_cath_annots.tsv" \
	--pfam-domains "$CdnPfamDir/domains.tsv" \
	--pfam-names "$CdnPfamDir/pfam_ref.tbl" \
	--motifs "$PalmDir/data/cdn_motifs.tsv" \
	--fasta "$PalmDir/data/cdn.fa" \
	--motif-style iii \
	--order "$PalmDir/data/cdn_figure_order.txt" \
	--scop-names "$PalmDir/data/scop_superfamily_names.txt" \
	--cath-names "$PalmDir/data/cath-b-newest-names" \
	--view-lo 1 \
	--output "$OutDir/cdn_cath_sf2.pdf"

echo "done: $OutDir/{rdrp,cdn}_cath_sf2.pdf"
