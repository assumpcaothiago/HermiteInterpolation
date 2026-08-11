PYTHON  ?= python3
LATEXMK ?= latexmk

DOCUMENT  := hermite_interpolation
TEX       := $(DOCUMENT).tex
PDF       := $(DOCUMENT).pdf
LOG       := $(DOCUMENT).log
DERIVATION := derive_hermite.py
GENERATED := generated_formulas.tex

LATEXMK_FLAGS := -pdf -interaction=nonstopmode -halt-on-error -file-line-error

.PHONY: all pdf formulas check clean distclean help

all: pdf

pdf: $(PDF)

$(PDF): $(TEX) $(GENERATED)
	$(LATEXMK) $(LATEXMK_FLAGS) $(TEX)

$(GENERATED): $(DERIVATION)
	$(PYTHON) $(DERIVATION) --emit-tex $(GENERATED)

formulas:
	$(PYTHON) $(DERIVATION) --emit-tex $(GENERATED)

check:
	$(PYTHON) $(DERIVATION) --check
	$(PYTHON) $(DERIVATION) --check-tex $(GENERATED)
	$(LATEXMK) $(LATEXMK_FLAGS) $(TEX)
	@! grep -E "(There were undefined references|Citation .* undefined|Reference .* undefined)" $(LOG)

# Keep the deliverable PDF and generated formula appendix.
clean:
	$(LATEXMK) -c $(TEX)

# Remove every reproducible artifact, including the delivered PDF and appendix.
distclean:
	$(LATEXMK) -C $(TEX)
	rm -f $(GENERATED)

help:
	@echo "make          regenerate formulas as needed and build $(PDF)"
	@echo "make formulas regenerate $(GENERATED)"
	@echo "make check    run symbolic checks, verify formulas, and build the PDF"
	@echo "make clean    remove auxiliary LaTeX files only"
	@echo "make distclean remove all reproducible artifacts"
