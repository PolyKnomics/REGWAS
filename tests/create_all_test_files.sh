#!/bin/bash
set -euo pipefail
set -x
cd "$(dirname "$0")"

pwd

curl  https://enkre.net/cgi-bin/code/bgen/raw/3ec770a829a753282b5cb45afc3f4eda036b246705b76f9037b6cc98c41a4194?at=example.8bits.bgen -o example.8bits.bgen
chmod -x example.8bits.bgen
curl https://enkre.net/cgi-bin/code/bgen/raw/3ec770a829a753282b5cb45afc3f4eda036b246705b76f9037b6cc98c41a4194?at=example.8bits.bgen.bgi -o example.8bits.bgen.bgi

R --vanilla --slave --no-save CMD BATCH create_small_phenotype.R 
R --vanilla --slave --no-save CMD BATCH create_big_phenotype.R
# R should not create these but remove them to be sure
rm -rf .RData
#put example files in right directory for example usage
mkdir -p ../inst/extdata/.
cp example.8bits* ../inst/extdata/.
cp *.tsv ../inst/extdata/.
