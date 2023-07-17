CXX = mpicxx

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = -g $(WARNINGS) -std=c++17 $(OPT) -fopenmp

LDLIBS = -L/sw/gcc/5.3.0/rtf/lib64 -L/home/joshua/usr/local/include/boost_lib/ -lboost_serialization -lboost_iostreams -lz

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3

COMMUTE_SRCS=Momentum.cpp Coefficient.cpp Operator.cpp Term.cpp WickTerm.cpp WickCleaner.cpp

HELPER_SRCS=PhaseHelper.cpp ModeHelper.cpp XPModes.cpp GeneralBasis.cpp
SQUARE_SRCS=HubbardCDW.cpp UsingBroyden.cpp SquareTripletPairing.cpp
CHAIN_SRCS=ChainTripletPairing.cpp
DOS_SRCS=BaseDOS.cpp Square.cpp SimpleCubic.cpp
SELFCON_SRCS=Selfconsistency/BroydenSolver.cpp
HBBRD_SRCS=$(addprefix Helper/, $(HELPER_SRCS)) $(SELFCON_SRCS) $(addprefix SquareLattice/, $(SQUARE_SRCS)) $(addprefix ChainLattice/, $(CHAIN_SRCS)) $(addprefix DensityOfStates/, $(DOS_SRCS)) ModelParameters.cpp
UTIL_SRCS=InputFileReader.cpp

PART_SRCS=Hubbard_Mean_Field.cpp
SRCS=$(addprefix Utility/, $(UTIL_SRCS)) $(addprefix Hubbard/, $(HBBRD_SRCS)) $(addprefix SymbolicOperators/, $(COMMUTE_SRCS)) $(PART_SRCS)

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: build build/main 

build/main: $(OBJS) | build
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: sources/%.cpp# sources/%.hpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

build:
	mkdir -p build
	mkdir -p build/Hubbard
	mkdir -p build/Hubbard/Helper
	mkdir -p build/Hubbard/Selfconsistency
	mkdir -p build/Hubbard/SquareLattice
	mkdir -p build/Hubbard/ChainLattice
	mkdir -p build/Hubbard/DensityOfStates
	mkdir -p build/Utility
	mkdir -p ../FermionCommute/build
	ln -s ../../FermionCommute/build build/SymbolicOperators

clean:
	rm -f $(OBJS)
	rm -f main
	rm -rf build
