#!/bin/bash -e

./hits_to_edf_and_topcat.bash lol_prefilter scop40 score
./hits_to_edf_and_topcat.bash lol_prefilter scop40x score

./hits_to_edf_and_topcat.bash lol_exhaustive scop40 score
./hits_to_edf_and_topcat.bash lol_exhaustive scop40x score
