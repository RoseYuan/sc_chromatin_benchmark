#!/bin/bash
sleep 5h
snakemake --configfile configs/Chen_2019.yaml --cores 1 --resources mem_mb=30000 --keep-going --rerun-triggers mtime |& tee -a std_output/output.Chen_2019.2.txt
