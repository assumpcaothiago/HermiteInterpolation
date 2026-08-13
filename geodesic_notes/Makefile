LATEXMK ?= latexmk

DOCUMENT := geodesics_in_adm_variables
TEX := $(DOCUMENT).tex
PDF := $(DOCUMENT).pdf
LATEXMK_FLAGS := -pdf -interaction=nonstopmode -halt-on-error -file-line-error

.PHONY: all pdf check clean distclean help

all: pdf

pdf: $(PDF)

$(PDF): $(TEX)
	$(LATEXMK) $(LATEXMK_FLAGS) $(TEX)

check: $(TEX)
	$(LATEXMK) $(LATEXMK_FLAGS) $(TEX)
	@! grep -E "(There were undefined references|Citation .* undefined|Reference .* undefined)" $(DOCUMENT).log

clean:
	$(LATEXMK) -c $(TEX)

distclean:
	$(LATEXMK) -C $(TEX)

help:
	@echo "make           build $(PDF)"
	@echo "make check     rebuild and reject unresolved references"
	@echo "make clean     remove auxiliary LaTeX files"
	@echo "make distclean also remove the generated PDF"
