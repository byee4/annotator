#!/bin/bash

gene_name2id.py \
--input inputs/shsap-holistic1.txt \
--gffdb inputs/gencode.v19.annotation.gtf.db \
--name_col Gene \
--output outputs/shsap-holistic1.geneids.txt
