#!/bin/bash -e

rm -f ../cache/reseek3*
rm -f ../edf/reseek3*
rm -f ../topcat/reseek3*

./cache_hits_reseek3.bash
./hits_to_edf_and_topcat_reseek3.bash
./summary.bash