#  Makefile for LaTeX documents
#
TARGET = shiva.ps
TOUCHFILE = shiva.tex
TOUCH=touch
PS = dvips
LATEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex
PS2PDF= ps2pdf
DELETE =  *.aux *.log *.ps *.dvi *.bbl *.blg *~ */*~ *.pdf
DVIPSFLAGS= -Pcmz -Pamz -Ppdf -G0 -tletter
TOUCH=touch
MAKE=make

.SUFFIXES: .pdf .tex .ps

all: $(TARGET)
	$(TOUCH) $(TOUCHFILE)
	$(MAKE) $(TARGET)

$(TARGET): 

.ps.pdf:
	$(PS2PDF) $<
	
.dvi.ps:
	$(PS) $(DVIPSFLAGS) -t letter -o $@ $<

.tex.dvi: 
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<

clean: 
	rm -f $(DELETE)
