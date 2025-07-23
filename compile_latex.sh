#!/bin/bash
CMD="latexmk -xelatex main.tex"
CMD+=" && latexmk -c main.tex"
CMD+=" && echo 'cleaning...' "
CMD+=" && rm -f *.bbl && rm *.xml"
eval $CMD