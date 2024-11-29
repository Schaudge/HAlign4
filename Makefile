###############################################################################
# Definitions
###############################################################################
FOLDER_WFA=PairwiseAlignment/WFA2-lib
FOLDER_LIB=PairwiseAlignment/WFA2-lib/lib

###############################################################################
# Flags & Folders
###############################################################################
FOLDER_BUILD=PairwiseAlignment/WFA2-lib/build
FOLDER_BUILD_CPP=PairwiseAlignment/WFA2-lib/build/cpp

UNAME=$(shell uname)

CC:=$(CC)
CPP:=$(CXX)

CC_FLAGS=-w -Wall -g -fPIE

AR=ar
AR_FLAGS=-rsc

GPP9_EXISTS := $(shell command -v g++-9 2> /dev/null)

ifeq ($(strip $(GPP9_EXISTS)),)
  CXX = g++
else
  CXX = g++-9
endif

###############################################################################
# Configuration rules
###############################################################################
LIB_WFA=$(FOLDER_LIB)/libwfa.a
LIB_WFA_CPP=$(FOLDER_LIB)/libwfacpp.a
SUBDIRS=PairwiseAlignment/WFA2-lib/alignment \
        PairwiseAlignment/WFA2-lib/bindings/cpp \
        PairwiseAlignment/WFA2-lib/system \
        PairwiseAlignment/WFA2-lib/utils \
        PairwiseAlignment/WFA2-lib/wavefront

all: CC_FLAGS+=-O3 -march=native #-flto -ffat-lto-objects
all: build h4

debug: build

ASAN_OPT=-fsanitize=address -fsanitize=undefined -fsanitize=shift -fsanitize=alignment
ASAN_OPT+=-fsanitize=signed-integer-overflow -fsanitize=bool -fsanitize=enum
ASAN_OPT+=-fsanitize=pointer-compare -fsanitize=pointer-overflow -fsanitize=builtin

# ASAN: ASAN_OPTIONS=detect_leaks=1:symbolize=1 LSAN_OPTIONS=verbosity=2:log_threads=1
asan: CC_FLAGS+=$(ASAN_OPT) -fno-omit-frame-pointer -fno-common
asan: build

###############################################################################
# Build rules
###############################################################################
build: setup
build: $(SUBDIRS) 
build: lib_wfa 

setup:
	@mkdir -p $(FOLDER_BUILD) $(FOLDER_BUILD_CPP) $(FOLDER_LIB)
    
lib_wfa: $(SUBDIRS)
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null
	$(AR) $(AR_FLAGS) $(LIB_WFA_CPP) $(FOLDER_BUILD)/*.o $(FOLDER_BUILD_CPP)/*.o 2> /dev/null

###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
    
.PHONY: $(SUBDIRS)

###############################################################################
# Rules
###############################################################################
LIBS=-fopenmp -lm
ifeq ($(UNAME), Linux)
  LIBS+=-lrt 
endif
        
h4: *.cpp $(LIB_WFA)
	g++ $(CC_FLAGS) -L$(FOLDER_LIB) -I$(FOLDER_WFA) \
	./PairwiseAlignment/NeedlemanWunshReusable.cpp \
	./SuffixArray/parallel_import.cpp \
	./Utils/Arguments.cpp \
	./Utils/Fasta.cpp \
	./Utils/Graph.cpp \
	./Utils/Insertion.cpp \
	./Utils/NucleicAcidColumn.cpp \
	./Utils/Utils.cpp \
	./multi-thread/multi.cpp \
	./StarAlignment/StarAligner.cpp \
	stmsa.cpp -o halign4 -static-libstdc++ -std=c++17 -lpthread -lwfacpp $(LIBS)

clean: 
	rm -rf $(FOLDER_BUILD) $(FOLDER_LIB) 2> /dev/null
	rm -rf $(FOLDER_TESTS)/*.alg $(FOLDER_TESTS)/*.log* 2> /dev/null
	rm h4
