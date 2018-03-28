#!/bin/bash

conda create -y -n annotator \
python=2.7 \
gffutils=0.8.7.1 \
bedtools=2.26 \
pybedtools=0.7.10 \
tqdm=4.14 \
future=0.16 \
futures=3.2.0 \
pytest=3.5.0

source activate annotator

conda install --override-channels \
-c conda-forge bzip2 # fixes some weird c-library issue

# already in install directory
python setup.py install

# adds cwl directory to path
export PATH=${PWD}/cwl:$PATH

