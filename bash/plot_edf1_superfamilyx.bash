#!/bin/bash -e

mkdir -p ../plots

inputs=()
shopt -s nullglob
for f in ../edf1/*.superfamilyx; do
	inputs+=(--input "$f")
done
if [ ${#inputs[@]} -eq 0 ]; then
	echo "no inputs match ../edf1/*.superfamilyx" >&2
	exit 1
fi

python ../py/plot_curve_edf1.py \
	"${inputs[@]}" \
	--colors ../info/colors.txt \
	--output ../plots/edf1_superfamilyx.svg
