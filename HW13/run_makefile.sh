#!/bin/bash

set -x
set -e

sample_id=$1
accession=$2

# Run the Makefile with the sample-specific parameters
make downloadsrr SAMPLE_ID=$sample_id SRR=$accession