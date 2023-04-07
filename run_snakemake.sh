#!/bin/bash
sleep 4h
snakemake --configfile configs/Ma_2019.yaml --cores 1 --resources mem_mb=30000 --keep-going --rerun-triggers mtime |& tee -a std_output/output.4.5.10.txt
