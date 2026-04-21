#!/bin/bash

source venv/bin/activate

smr="data/smr"
myeqtl="data/westra_eqtl_hg19"
myeqtlbig="data/eqtlgen-sparse"

# Create the eqtlgen dataset
# Note, takes a while
wget https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz
tar xzvf cis-eQTL-SMR_20191212.tar.gz
gunzip cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense*
${smr} --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --make-besd --out ${myeqtlbig}

# Use SMR
time $smr \
    --beqtl-summary $myeqtl \
    --query 5e-5 \
    --snp-chr 1 \
    --from-snp-kb 100 \
    --to-snp-kb 2000 \
    --probe-chr 1 \
    --from-probe-kb 1000 \
    --to-probe-kb 2000 \
    --out results/myquery_smr

# This takes too long to run?!
time $smr \
    --beqtl-summary $myeqtlbig \
    --query 5e-5 \
    --snp-chr 1 \
    --from-snp-kb 100 \
    --to-snp-kb 2000 \
    --probe-chr 1 \
    --from-probe-kb 1000 \
    --to-probe-kb 2000 \
    --out results/myquery_smr

# Use BESD format with besdq
source venv/bin/activate

time python3 -m besdq.cli \
    --beqtl-summary $myeqtl \
    --query 5.0e-4 \
    --snp-chrpos 1:1191870 \
    --probe-chrpos 1:1140818 \
    --out results/remyquery

pip install -e .

# Use installed package
time besdq \
    --beqtl-summary $myeqtl \
    --query 5e-4 \
    --snp-chr 1 \
    --from-snp-kb 100 \
    --to-snp-kb 2000 \
    --probe-chr 1 \
    --from-probe-kb 1000 \
    --to-probe-kb 2000 \
    --out results/myquery

time besdq \
    --beqtl-summary $myeqtlbig \
    --query 5e-4 \
    --snp-chr 1 \
    --from-snp-kb 100 \
    --to-snp-kb 2000 \
    --probe-chr 1 \
    --from-probe-kb 1000 \
    --to-probe-kb 2000 \
    --out results/myquery

# Use BESD-index format with besdq
besdq --beqtl-summary data/westra_eqtl_hg19 --index data/westra_eqtl_hg19.db

time besdq \
    --besd-index ${myeqtl}.db \
    --query 5e-4 \
    --snp-chr 1 \
    --from-snp-kb 100 \
    --to-snp-kb 2000 \
    --probe-chr 1 \
    --from-probe-kb 1000 \
    --to-probe-kb 2000 \
    --out results/myquery

time besdq \
    --besd-index ${myeqtlbig}.db \
    --query 5e-4 \
    --snp-chr 1 \
    --from-snp-kb 100 \
    --to-snp-kb 2000 \
    --probe-chr 1 \
    --from-probe-kb 1000 \
    --to-probe-kb 2000 \
    --out results/myquery
