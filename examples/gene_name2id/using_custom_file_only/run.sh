#!/bin/bash

gene_name2id.py \
--input inputs/shsap-holistic1.txt \
--custom inputs/ensembl2entrez_GRCh37_mart_export.tsv \
--name_col EntrezID \
--output outputs/shsap-holistic1.geneids.txt
