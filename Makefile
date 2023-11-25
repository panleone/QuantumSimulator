CC=gcc
CXX=g++
CFLAGS=-O3 -march=native -ffast-math
CXXFLAGS=$(CFLAGS)
LDFLAGS=-llapack -lopenblas -lfftw3_threads -lfftw3 -lm
SRC=$(shell find . -path ./unitTest -prune -o -name '*.cpp' ! -name 'main.cpp' -print)
HEADERS = $(shell find . -type f -name '*.h')
OBJ = $(SRC:%.cpp=%.o)
BIN = quantumSimulator

$(BIN): $(OBJ) main.o
	$(CXX) $(CFLAGS) $(OBJ) main.o $(LIBS) $(LDFLAGS) -o $@

%.o: %.c Makefile
	$(CC) $(CFLAGS) $? -c -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $? -c -o $@

.PHONY: clean format test
clean:
	find . -type f -name '*.o*' -exec rm {} +
	-rm "$(BIN)"

format:
	clang-format -i $(SRC) $(HEADERS)
SRC_TEST=$(shell find unitTest -iname '*.cpp')
OBJ_TEST=$(SRC_TEST:%.cpp=%.o)
.PHONY: test
test: $(OBJ) $(OBJ_TEST)
	$(CXX) $(CFLAGS) $(OBJ) $(OBJ_TEST) $(LIBS) $(LDFLAGS) -o $@
