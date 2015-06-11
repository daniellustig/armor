#!/bin/bash

set -e

mkdir -p outputs

echo "Generating shim FSMs..."
./report.py -o outputs -l outputs/log.txt

echo "Generating summary reports..."
for i in isca-2x2 isca-2x4 isca-4x4 isca-4x6 normal
do
  chmod +x ./outputs/$i/build_graphs.sh
  ./outputs/$i/build_graphs.sh
  cd outputs/$i && pdflatex summary.tex && cd ../..
done

echo "Done!"
echo
echo
for i in isca-2x2 isca-2x4 isca-4x4 isca-4x6 normal
do
  echo Report available at outputs/$i/summary.pdf
done

echo "Building tech report..."
pdftk doc/dlustig_ISCA15.pdf outputs/isca-4x4/summary.pdf cat output doc/dlustig_ISCA15_TR.pdf
echo "Tech report available at doc/dlustig_ISCA15_TR.pdf"
