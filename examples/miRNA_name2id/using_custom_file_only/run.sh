#!/bin/bash

miRNA_name2id.py \
--input inputs/all_mirnas.csv \
--sep , \
--custom inputs/ensembl2name_GCm38_mart_export.tsv \
--name_col miRNA \
--output outputs/all_mirnas.with_ids.csv
