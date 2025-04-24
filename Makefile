all: main
clean:
	-@rm -f main && rm -rf bin/

CXX = g++ -std=c++20 -W -Wall -march=native -Ifsst -Isrc/ -I.
OPT=-O3 -fno-stack-protector
ifdef NDEBUG
	OPT += -DNDEBUG
endif

COMPILER_VERSION = $(shell ${CXX} --version)
ifneq ($(findstring clang,$(COMPILER_VERSION)),)
   CXX += -Wno-unqualified-std-cast-call
endif

# Build/link the main.
main: main.cpp fsst/build/libfsst.a bin/BenchmarkDriver.o bin/FsstWrapper.o bin/StateMachine.o src/algos/StdFind.hpp src/algos/StartsWith.hpp src/algos/Comet.hpp
	$(CXX) -o main $(OPT) main.cpp -flto bin/BenchmarkDriver.o bin/FsstWrapper.o bin/StateMachine.o -Lfsst/build -lfsst 

# Build for individual files.
bin/BenchmarkDriver.o: src/BenchmarkDriver.cpp src/BenchmarkDriver.hpp
	@mkdir -p bin
	$(CXX) -c $(OPT) src/BenchmarkDriver.cpp -o bin/BenchmarkDriver.o

bin/FsstWrapper.o: src/FsstWrapper.cpp src/FsstWrapper.hpp
	@mkdir -p bin
	$(CXX) -c $(OPT) src/FsstWrapper.cpp -o bin/FsstWrapper.o

bin/StateMachine.o: src/StateMachine.hpp src/StateMachine.cpp
	@mkdir -p bin
	$(CXX) -c $(OPT) src/StateMachine.cpp -o bin/StateMachine.o