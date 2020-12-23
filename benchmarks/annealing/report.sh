#!/bin/bash

python $PWD/plot.py
pdflatex $PWD/report.tex
pdflatex $PWD/report.tex

rm -f $PWD/HYP-24_chimera_*.pdf
rm -f *.aux *.log *.out
