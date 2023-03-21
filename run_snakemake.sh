#!/bin/bash

snakemake --configfile configs/candidate1.yaml --cores 12 --resources mem_mb=30000 |& tee -a output.txt
