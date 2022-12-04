This work is licensed under a Creative Commons ‘Attribution-NonCommercial-ShareAlike 4.0 International’ licence. 
More information in LICENSE.txt.

To compile the LaTeX into a PDF I recommend a full install of texlive (texlive-full) on Linux. 
If you use Mac or Windows products: ... well there's your problem, you shouldn't! You have to figure it out yourself!

For a complete compile use this command sequence in your terminal while in the base folder of the project:
pdflatex -synctex=1 -interaction=nonstopmode thesis
pdflatex -synctex=1 -interaction=nonstopmode thesis
biber thesis
pdflatex -synctex=1 -interaction=nonstopmode thesis
makeglossaries thesis
pdflatex -synctex=1 -interaction=nonstopmode thesis
biber thesis
pdflatex -synctex=1 -interaction=nonstopmode thesis
pdflatex -synctex=1 -interaction=nonstopmode thesis

If you change Feynman diagrams you need to compile them again with:
cd feynman
mpost *.mp
cd ..
pdflatex -synctex=1 -interaction=nonstopmode thesis

Some graphs are drawn using ROOT. To redraw them you need a full ROOT install (especially the OpenGL parts for colour alpha). To install root go to https://root.cern.ch/.
The scripts can then be executed by typing the following into your terminal
cd code 
root -l NameOfScript.C

In order to adapt other graphics, edit the .svg files in Inkscape and thereafter convert the files to PDF format (only the PDFs are included in the thesis). 
