#!/usr/bin/env bash

# Project name
PROJ=$1

# Input directory with "codes" files
IDIR="intermediate/projects/${PROJ}"

# Results output directory for colsum and poolcount files
ODIR="results/projects/${PROJ}"

# Path to FEBA combineBarSeq.pl script
CBS="/ssd/common/tools/feba/bin/combineBarSeq.pl"

# Run combineBarSeq.pl script on "codes" files
${CBS} ${ODIR}/${PROJ} data/projects/${PROJ}.poolfile.tab ${IDIR}/*.codes
