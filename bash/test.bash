#!/bin/bash -e

cd ../test

python ../py/build_derived_info.py \
	--lookup test.lookup \
	--output derived_info.json
