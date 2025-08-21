# Customize compiler version here if desired:
CC = gcc
CXX = g++

#
RM = rm -f
LN = ln -f

CXXFLAGS+= -fopenmp
CXXFLAGS+= -Iextern
CXXFLAGS+= -std=c++17
CXXFLAGS+= -ffast-math -funsafe-math-optimizations -fno-finite-math-only -O3
CXXFLAGS+= -DNDEBUG
CXXFLAGS+= -Wno-unused-variable -Wno-maybe-uninitialized

MAIN=amplitrop

.PHONY: depend clean

all:    $(MAIN)

.depend :

amplitrop : amplitrop.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

amplitrop.o : amplitrop.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


depend: .depend

.depend : amplitrop.cpp
	$(RM) ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(MAIN) amplitrop.o .depend

include .depend

