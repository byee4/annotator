#!/bin/bash

miRNA_name2id.py \
--input inputs/all_mirnas.csv \
--sep , \
--gffdb inputs/mmu.mirBase_v21.GCm38.gff3.db \
--custom inputs/custom_name2acc_mirBase_v21.GCm38.tsv \
--name_col miRNA \
--output outputs/all_mirnas.with_ids.csv \
--add_mature
