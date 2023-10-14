# By default, this Makefile only builds the `fastdna` binary. If you want to
# build the paper's PDF document (and its prerequisite experimental data),
# use `make cal.pdf` (or whatever document name you want
#
CXX=g++-13
OPTFLAGS=-O0 -g
OPTFLAGS=-O3 -DNDEBUG
CXXFLAGS=-std=c++2a -Wall -Wextra -pedantic -mbmi2 -mtune=native -pthread

OUTPUTS=*.pdf *.tex _main.* libs _book _bookdown_files *.html *.log *.aux *.bbl *.blg *.dvi *.out *.xcp
R=Rscript

.PHONY: all

all: fastdna

fastdna: fastdna.o
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) -o $@ $^ -lrt -flto

%.o: %.cc
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) -c  $^

%.pdf: %.Rmd %.csv Makefile fastdna.bib
	if [ -f _main.Rmd ]; then \
		rm _main.Rmd; \
	fi
	Rscript -e 'library(bookdown); bookdown::render_book("$(basename $<).Rmd", "bookdown::pdf_book", config_file = "$(basename $<).yml")'
	mv _book/_main.pdf $@


%.csv: fastdna
	python3 $(basename $@)-experiments.py


clean:
	rm -rf fastdna *.o *.dSYM _main.*
	$(R) -e 'bookdown::clean_book(clean = T)'
	rm -rf $(OUTPUTS)
