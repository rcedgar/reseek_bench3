#!/bin/bash -e

./stage1a_fast_hits_to_edf.bash
./stage1b_fast_hits_to_topcat.bash

./summary.bash

echo SECONDS=$SECONDS
