#!/bin/bash

snakemake --configfile configs/Chen_2019.yaml --cores 3 --resources mem_mb=30000 --keep-going=true |& tee -a output.txt
