CXX = mpicxx

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = $(WARNINGS) -std=c++20 $(OPT) -fopenmp

LDL_DIRS = -L/sw/gcc/5.3.0/rtf/lib64 -L/home/joshua/usr/local/include/boost_lib/
LDLIBS = $(LDL_DIRS) -lboost_serialization -lboost_iostreams -lz ~/usr/local/lib/libSymbolicOperators.a

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3# -ffast-math

HELPER_SRCS=PhaseHelper.cpp Plaquette.cpp ModeHelper.cpp XPModes.cpp GeneralBasis.cpp TermOnSquare.cpp SquareXP.cpp SquareGeneral.cpp
SQUARE_SRCS=HubbardCDW.cpp UsingBroyden.cpp SquareTripletPairing.cpp
CHAIN_SRCS=ChainTripletPairing.cpp
DOS_SRCS=BaseDOS.cpp Square.cpp SimpleCubic.cpp
MODEL_SRCS=$(addprefix SquareLattice/, $(SQUARE_SRCS)) $(addprefix ChainLattice/, $(CHAIN_SRCS)) ModelParameters.cpp EMCoupling.cpp
HBBRD_SRCS=$(addprefix Helper/, $(HELPER_SRCS)) $(addprefix DensityOfStates/, $(DOS_SRCS)) $(addprefix Models/, $(MODEL_SRCS))

PART_SRCS=Handler/HandlerBase.cpp Handler/TestHandler.cpp Handler/ModeHandler.cpp Handler/PhaseHandler.cpp Handler/UnknownBoundaryHandler.cpp Handler/ModeDispersionHandler.cpp HubbardMeanField.cpp
SRCS=$(addprefix Hubbard/, $(HBBRD_SRCS)) $(PART_SRCS)

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: build build/main 

debug: CXXFLAGS += -g -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment
debug: build build/main

build/main: $(OBJS) | build 
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: sources/%.cpp# sources/%.hpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

build:
	mkdir -p build
	mkdir -p build/Handler
	mkdir -p build/Hubbard
	mkdir -p build/Hubbard/Helper
	mkdir -p build/Hubbard/Models
	mkdir -p build/Hubbard/Models/SquareLattice
	mkdir -p build/Hubbard/Models/ChainLattice
	mkdir -p build/Hubbard/DensityOfStates

clean:
	rm -rf build
