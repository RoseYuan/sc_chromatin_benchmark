#!/bin/bash

snakemake --configfile configs/candidate1.yaml --cores 3 --resources mem_mb=30000 --keep-going True |& tee -a output.txt
