# Makefile for homval C++ tools (build under WSL/Linux).
#
#   make            # build hits_to_edf, hits_to_bin, hits_to_topcat
#   make test       # compare against golden .edf and .tcat fixtures
#   make clean

CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -pthread -Wall -Wextra -pedantic
LDFLAGS ?= -pthread

BINDIR := bin
EDF := $(BINDIR)/hits_to_edf
BIN := $(BINDIR)/hits_to_bin
TOPCAT := $(BINDIR)/hits_to_topcat

.PHONY: all test clean

all: $(EDF) $(BIN) $(TOPCAT)

$(BINDIR):
	mkdir -p $(BINDIR)

$(EDF): hits_to_edf.cpp hits_bin.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) hits_to_edf.cpp -o $(EDF) $(LDFLAGS)

$(BIN): hits_to_bin.cpp hits_bin.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) hits_to_bin.cpp -o $(BIN) $(LDFLAGS)

$(TOPCAT): hits_to_topcat.cpp hits_bin.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) hits_to_topcat.cpp -o $(TOPCAT) $(LDFLAGS)

# Regression: text and bin pipelines must match test/expected.* fixtures
test: $(EDF) $(BIN) $(TOPCAT)
	@set -e; \
	for truth in superfamily superfamilyx; do \
	  out=$$(mktemp); \
	  ./$(EDF) \
	    --hits ../test/test_hits.tsv \
	    --lookup ../test/test.lookup \
	    --derived-info ../test/derived_info.json \
	    --truth $$truth \
	    --score \
	    --algo brute_force_reference \
	    --output $$out; \
	  if ! diff -u --strip-trailing-cr ../test/expected.$$truth.edf $$out; then \
	    echo "FAIL text edf $$truth"; rm -f $$out; exit 1; \
	  fi; \
	  echo "OK text edf $$truth"; \
	  rm -f $$out; \
	done; \
	for truth in topsf topfold; do \
	  out=$$(mktemp); \
	  ./$(TOPCAT) \
	    --hits ../test/test_hits.tsv \
	    --lookup ../test/test.lookup \
	    --derived-info ../test/derived_info.json \
	    --truth $$truth \
	    --score \
	    --algo brute_force_reference \
	    --output $$out; \
	  if ! diff -u --strip-trailing-cr ../test/expected.$$truth.tcat $$out; then \
	    echo "FAIL text topcat $$truth"; rm -f $$out; exit 1; \
	  fi; \
	  echo "OK text topcat $$truth"; \
	  rm -f $$out; \
	done; \
	cache=$$(mktemp); \
	./$(BIN) \
	  --hits ../test/test_hits.tsv \
	  --lookup ../test/test.lookup \
	  --output $$cache.bin; \
	for truth in superfamily superfamilyx; do \
	  out=$$(mktemp); \
	  ./$(EDF) \
	    --bin $$cache.bin \
	    --lookup ../test/test.lookup \
	    --derived-info ../test/derived_info.json \
	    --truth $$truth \
	    --score \
	    --algo brute_force_reference \
	    --output $$out; \
	  if ! diff -u --strip-trailing-cr ../test/expected.$$truth.edf $$out; then \
	    echo "FAIL bin edf $$truth"; rm -f $$out $$cache.bin; exit 1; \
	  fi; \
	  echo "OK bin edf $$truth"; \
	  rm -f $$out; \
	done; \
	for truth in topsf topfold; do \
	  out=$$(mktemp); \
	  ./$(TOPCAT) \
	    --bin $$cache.bin \
	    --lookup ../test/test.lookup \
	    --derived-info ../test/derived_info.json \
	    --truth $$truth \
	    --score \
	    --algo brute_force_reference \
	    --output $$out; \
	  if ! diff -u --strip-trailing-cr ../test/expected.$$truth.tcat $$out; then \
	    echo "FAIL bin topcat $$truth"; rm -f $$out $$cache.bin; exit 1; \
	  fi; \
	  echo "OK bin topcat $$truth"; \
	  rm -f $$out; \
	done; \
	prefix=$$(mktemp); \
	./$(EDF) \
	  --bin $$cache.bin \
	  --lookup ../test/test.lookup \
	  --derived-info ../test/derived_info.json \
	  --all-truths \
	  --score \
	  --algo brute_force_reference \
	  --output $$prefix; \
	for truth in family superfamily superfamilyx fold; do \
	  if [ "$$truth" = "superfamily" ] || [ "$$truth" = "superfamilyx" ]; then \
	    if ! diff -u --strip-trailing-cr ../test/expected.$$truth.edf $$prefix.$$truth; then \
	      echo "FAIL all-truths edf $$truth"; rm -f $$cache.bin $$prefix.*; exit 1; \
	    fi; \
	    echo "OK all-truths edf $$truth"; \
	  fi; \
	done; \
	tprefix=$$(mktemp); \
	./$(TOPCAT) \
	  --bin $$cache.bin \
	  --lookup ../test/test.lookup \
	  --derived-info ../test/derived_info.json \
	  --all-truths \
	  --score \
	  --algo brute_force_reference \
	  --output $$tprefix; \
	for truth in topsf topfold; do \
	  if ! diff -u --strip-trailing-cr ../test/expected.$$truth.tcat $$tprefix.$$truth; then \
	    echo "FAIL all-truths topcat $$truth"; rm -f $$cache.bin $$prefix.* $$tprefix.*; exit 1; \
	  fi; \
	  echo "OK all-truths topcat $$truth"; \
	done; \
	rm -f $$cache.bin $$prefix.* $$tprefix.*

clean:
	rm -f $(EDF) $(BIN) $(TOPCAT)
