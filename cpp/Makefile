# Makefile for homval C++ tools (build under WSL/Linux).
#
#   make            # build hits_to_edf, hits_to_bin, hits_to_topcat, build_derived_info
#   make test       # compare against golden .edf and .tcat fixtures
#   make clean

CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -pthread -Wall -Wextra -pedantic
LDFLAGS ?= -pthread

BINDIR := bin
EDF := $(BINDIR)/hits_to_edf
BIN := $(BINDIR)/hits_to_bin
TOPCAT := $(BINDIR)/hits_to_topcat
DERIVED := $(BINDIR)/build_derived_info

HEADERS := hits_bin.hpp taxonomy.hpp derived_info.hpp

.PHONY: all test clean

all: $(EDF) $(BIN) $(TOPCAT) $(DERIVED)

$(BINDIR):
	mkdir -p $(BINDIR)

$(EDF): hits_to_edf.cpp $(HEADERS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) hits_to_edf.cpp -o $(EDF) $(LDFLAGS)

$(BIN): hits_to_bin.cpp $(HEADERS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) hits_to_bin.cpp -o $(BIN) $(LDFLAGS)

$(TOPCAT): hits_to_topcat.cpp $(HEADERS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) hits_to_topcat.cpp -o $(TOPCAT) $(LDFLAGS)

$(DERIVED): build_derived_info.cpp $(HEADERS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) build_derived_info.cpp -o $(DERIVED) $(LDFLAGS)

# Regression: text and bin pipelines must match test/expected.* fixtures
test: $(EDF) $(BIN) $(TOPCAT) $(DERIVED)
	@set -e; \
	./$(DERIVED) --lookup ../test/test.lookup --output ../test/derived_info.json; \
	for truth in superfamily superfamilyx; do \
	  out=$$(mktemp); \
	  ./$(EDF) \
	    --hits ../test/test_hits.tsv \
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
	rm -f $(EDF) $(BIN) $(TOPCAT) $(DERIVED)
