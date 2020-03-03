#!/bin/bash

GMR=$(awk -v OFS='\t' '{sum+=$5}END{print sum}' ${id}.collapsed.bed)
MGMR=$(awk '{print +$1 }' <<< $GMR)
MGMR=$(bc -l <<< $MGMR/1000000)
echo $MGMR
awk -v OFS='\t' -v MGMR="$MGMR" '{print $1,$2,$3,$4,$5/MGMR,$6,$7,$8}' ${id}.collapsed.bed > ${id}.norm.bed



