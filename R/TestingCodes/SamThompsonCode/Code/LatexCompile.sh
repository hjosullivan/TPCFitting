#!/bin/bash
# Author: Sam Thompson samuel.thompson14@imperial.ac.uk
# Script: LatexCompile.sh
# Desc: Compiles latex file which it is directed to
# Arguments: input_latex_file_parent_directory, file name, output_file_directory
# Date: Nov 2014
cd "$1"
pdflatex "$2".tex
pdflatex "$2".tex
bibtex "$2"
bibtex "$2"
pdflatex "$2".tex
pdflatex "$2".tex

##Cleanup
rm *~
rm "$1/"*.aux
rm "$1/"*.dvi
rm "$1/"*.log
rm "$1/"*.nav
rm "$1/"*.out
rm "$1/"*.snm
rm "$1/"*.toc
rm "$1/"*.bbl
rm "$1/"*.blg
mv "$1/$2.pdf" "$3"

exit

