# Makefile to compile all .tex file in the current directory (and recursively)
# with latexmk. Supports parallel compilation and TEX compiler directive
#
# Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
# Created: Sun 1 Dec 2019, flying to Seattle

LC    = latexmk
LOGS  = .auxtex
FLAGS = -quiet -outdir=$(LOGS) -M -MP -MF $(LOGS)/$*.deps
SRCS  = $(wildcard *.tex)
PDFS  = $(patsubst %.tex,%.pdf,$(SRCS))

all : $(PDFS)

%.pdf : %.tex $(LOGS) FORCE_MAKE
	@if [[ $$(head -n1 $*.tex | grep -oh "\w*tex\w*") == 'xelatex' ]]; then \
	    echo [$(shell pwd)] $(LC) $(FLAGS) -pdfxe -pdfxelatex=\"xelatex -halt-on-error\" $<; \
	    eval $(LC) $(FLAGS) -pdfxe -pdfxelatex=\"xelatex -halt-on-error\" $< &> /dev/null; \
	    exit $$?; \
	elif [[ $$(head -n1 $*.tex | grep -oh "\w*tex\w*") == 'lualatex' ]]; then \
	    echo [$(shell pwd)] $(LC) $(FLAGS) -pdflua -pdflualatex=\"lualatex -halt-on-error\" $<; \
	    eval $(LC) $(FLAGS) -pdflua -pdflualatex=\"lualatex --halt-on-error\" $< &> /dev/null; \
	    exit $$?; \
	else \
	    echo [$(shell pwd)] $(LC) $(FLAGS) -pdf -pdflatex=\"pdflatex -halt-on-error\" $<; \
	    eval $(LC) $(FLAGS) -pdf -pdflatex=\'pdflatex -halt-on-error\' $< &> /dev/null; \
	    exit $$?; \
	fi
	@cp $(LOGS)/$@ .

$(LOGS) :
	@echo [$(shell pwd)] mkdir -p $(LOGS)
	@mkdir -p $(LOGS)

clean :
	@echo [$(shell pwd)] rm -rf $(PDFS) $(LOGS)
	-@rm -rf $(PDFS) $(LOGS)

cleanall : clean
	+cd img && make clean; cd ..

.PHONY : FORCE_MAKE clean
-include $(LOGS)/*.deps
