# By default, this Makefile only builds the `fastdna` binary. Use `make cal.pdf`
# to build the paper from the checked-in experiment data, or `make cal-data` to
# rerun the experiments and regenerate that data.
#
CXX=g++-13
OPTFLAGS=-O0 -g
OPTFLAGS=-O3 -DNDEBUG
CXXFLAGS=-std=c++2a -Wall -Wextra -pedantic -mbmi2 -mtune=native -pthread

OUTPUTS=*.pdf main.tex _main.* libs _book _bookdown_files *.html *.log *.aux *.bbl *.blg *.dvi *.out *.xcp
UVR_VERSION=v0.4.7
UVR_DIR=$(HOME)/.cache/fastdna/uvr/$(UVR_VERSION)
UVR=$(UVR_DIR)/uvr
UVR_LIBRARY=$(HOME)/.cache/fastdna/$(notdir $(CURDIR))/r-library

.PHONY: all r-deps cal-data clean

all: fastdna cal.pdf

$(UVR):
	mkdir -p $(UVR_DIR)
	curl -fsSL https://raw.githubusercontent.com/nbafrank/uvr/$(UVR_VERSION)/install.sh -o $(UVR_DIR)/install-uvr.sh
	UVR_VERSION=$(UVR_VERSION) UVR_INSTALL_DIR=$(abspath $(UVR_DIR)) sh $(UVR_DIR)/install-uvr.sh
	rm -f $(UVR_DIR)/install-uvr.sh

r-deps: $(UVR) uvr.toml uvr.lock
	UVR_LIBRARY=$(UVR_LIBRARY) $(UVR) sync --frozen --unattended --library $(UVR_LIBRARY)

fastdna: fastdna.o
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) -o $@ $^ -lrt -flto

%.o: %.cc
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) -c  $^

%.pdf: %.Rmd %.csv Makefile fastdna.bib %.yml cal-build.R uvr.toml uvr.lock | r-deps
	UVR_LIBRARY=$(UVR_LIBRARY) $(UVR) run --no-companion --unattended cal-build.R -- $(basename $<)
	mv _book/_main.pdf $@


clean: | r-deps
	rm -rf fastdna *.o *.dSYM _main.*
	UVR_LIBRARY=$(UVR_LIBRARY) $(UVR) run --no-companion --unattended cal-build.R -- clean
	rm -rf $(OUTPUTS)

cal-data: fastdna
	python3 cal-experiments.py
