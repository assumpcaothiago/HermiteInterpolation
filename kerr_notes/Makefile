PYTHON  ?= python3
LATEXMK ?= latexmk

DOCUMENT := kerr_quasi_isotropic_note
TEX      := $(DOCUMENT).tex
PDF      := $(DOCUMENT).pdf
LOG      := $(DOCUMENT).log
CHECKER  := kerr_quasi_isotropic_check.py

LATEXMK_FLAGS := -pdf -interaction=nonstopmode -halt-on-error -file-line-error

.PHONY: all pdf check clean distclean help

all: pdf

pdf: $(PDF)

$(PDF): $(TEX)
	$(LATEXMK) $(LATEXMK_FLAGS) $(TEX)

check: $(TEX) $(CHECKER)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m doctest $(CHECKER)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) $(CHECKER)
	$(LATEXMK) $(LATEXMK_FLAGS) $(TEX)
	@! grep -Eq "undefined references|Citation .* undefined|Reference .* undefined" $(LOG)

clean:
	$(LATEXMK) -c $(TEX)
	$(RM) -r __pycache__

distclean:
	$(LATEXMK) -C $(TEX)
	$(RM) -r __pycache__

help:
	@echo "make           build $(PDF)"
	@echo "make check     run symbolic checks and build the PDF"
	@echo "make clean     remove auxiliary LaTeX files and Python bytecode"
	@echo "make distclean also remove the generated PDF"
