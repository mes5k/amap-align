################################################################################
# Makefile for amap
################################################################################

################################################################################
# 1) Choose C++ compiler.
################################################################################

CXX = g++

################################################################################
# 2) Set C++ flags.
#    a) DEBUG mode -- no optimizations, enable SafeVector checking, no inlining
#    b) PROFILE mode -- for gprof
#    c) RELEASE mode
################################################################################

OTHERFLAGS = -DNumInsertStates=1 -DVERSION='"AMAP.2.2"'

# debug mode    
#CXXFLAGS = -g -W -Wall -pedantic -DENABLE_CHECKS -fno-inline $(OTHERFLAGS)
#CXXFLAGS = -g -W -Wall -pedantic -fno-inline $(OTHERFLAGS)

# profile mode
#CXXFLAGS = -pg -W -Wall -pedantic $(OTHERFLAGS)

# release mode
#CXXFLAGS = -O3 -W -Wall -pedantic -DNDEBUG $(OTHERFLAGS) -mmmx -msse -msse2 -mfpmath=sse -march=pentium4 -mcpu=pentium4 -funroll-loops -fomit-frame-pointer 
CXXFLAGS = -O3 -W -Wall -pedantic -DNDEBUG $(OTHERFLAGS) -funroll-loops

################################################################################
# 3) Dependencies
################################################################################

TARGETS = amap

.PHONY : all
all : $(TARGETS)

amap : MultiSequenceDag.h MultiSequence.h ProbabilisticModel.h ScoreType.h Sequence.h FileBuffer.h SparseMatrix.h EvolutionaryTree.h Defaults.h SafeVector.h Amap.cc
	$(CXX) $(CXXFLAGS) -lm -o amap Amap.cc 

.PHONY : clean
clean:
	rm -f $(TARGETS)
