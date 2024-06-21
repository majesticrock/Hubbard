CXX = mpicxx

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include -I ../../../PhdUtility/Utility/include -I ../../../PhdUtility/FermionCommute/sources

CXXFLAGS = $(WARNINGS) -std=c++20 $(OPT) -fopenmp

LDLIBS = -L/sw/gcc/5.3.0/rtf/lib64 -L/home/joshua/usr/local/include/boost_lib/ -lboost_serialization -lboost_iostreams -lz

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3# -ffast-math

COMMUTE_SRCS=TermLoader.cpp Coefficient.cpp IndexWrapper.cpp Momentum.cpp MomentumList.cpp Operator.cpp Term.cpp OperatorType.cpp WickOperator.cpp WickOperatorTemplate.cpp WickTerm.cpp WickSymmetry.cpp Wick.cpp

HELPER_SRCS=PhaseHelper.cpp Plaquette.cpp ModeHelper.cpp XPModes.cpp GeneralBasis.cpp TermOnSquare.cpp SquareXP.cpp SquareGeneral.cpp
SQUARE_SRCS=HubbardCDW.cpp UsingBroyden.cpp SquareTripletPairing.cpp
CHAIN_SRCS=ChainTripletPairing.cpp
DOS_SRCS=BaseDOS.cpp Square.cpp SimpleCubic.cpp
HBBRD_SRCS=$(addprefix Helper/, $(HELPER_SRCS)) $(addprefix SquareLattice/, $(SQUARE_SRCS)) $(addprefix ChainLattice/, $(CHAIN_SRCS)) $(addprefix DensityOfStates/, $(DOS_SRCS)) ModelParameters.cpp EMCoupling.cpp

PART_SRCS=Handler/HandlerBase.cpp Handler/TestHandler.cpp Handler/ModeHandler.cpp Handler/PhaseHandler.cpp Handler/UnknownBoundaryHandler.cpp Handler/ModeDispersionHandler.cpp HubbardMeanField.cpp
SRCS=$(addprefix Hubbard/, $(HBBRD_SRCS)) $(addprefix SymbolicOperators/, $(COMMUTE_SRCS)) $(PART_SRCS)

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: sources/SymbolicOperators build build/main 

debug: CXXFLAGS += -g -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment
debug: build build/main

build/main: $(OBJS) | build 
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: sources/%.cpp# sources/%.hpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

sources/SymbolicOperators:
	ln -s ../../../../PhdUtility/FermionCommute/sources/SymbolicOperators sources/SymbolicOperators

build:
	mkdir -p build
	mkdir -p build/Handler
	mkdir -p build/Hubbard
	mkdir -p build/Hubbard/Helper
	mkdir -p build/Hubbard/SquareLattice
	mkdir -p build/Hubbard/ChainLattice
	mkdir -p build/Hubbard/DensityOfStates
	mkdir -p build/SymbolicOperators

clean:
	rm -rf build
