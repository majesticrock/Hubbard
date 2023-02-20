CXX = mpicxx

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = -g $(WARNINGS) -std=c++17 $(OPT) -fopenmp -DNDEBUG

LDLIBS = ~/usr/local/include/boost_lib/libboost_serialization.a

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3

COMMUTE_SRCS=Momentum.cpp Coefficient.cpp Operator.cpp Term.cpp WickTerm.cpp WickCleaner.cpp

HBBRD_SRCS=Model.cpp ModelSubClasses.cpp BasicHubbardModel.cpp HubbardCDW.cpp UsingBroyden.cpp
UTIL_SRCS=InputFileReader.cpp OutputWriter.cpp Roots_Broyden.cpp Resolvent.cpp

PART_SRCS=
SRCS=$(addprefix Utility/, $(UTIL_SRCS)) $(addprefix Hubbard/, $(HBBRD_SRCS)) $(addprefix SymbolicOperators/, $(COMMUTE_SRCS)) $(PART_SRCS) Hubbard_Mean_Field.cpp

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: build build/main 

build/main: $(OBJS) | build
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: sources/%.cpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

build:
	mkdir -p build
	mkdir -p build/Hubbard
	mkdir -p build/Utility
	mkdir -p build/SymbolicOperators

clean:
	rm -f $(OBJS)
	rm -f main
	rm -rf build
