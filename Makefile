SHELL := /bin/bash

OPT=-O3 -DNDEBUG

all: main
clean:
	-@rm -f main

main: main.cc fsst/build/libfsst.a
	g++ -std=c++20 -W -Wall -o main $(OPT) main.cc -llz4 -Ifsst -Lfsst/build -lfsst 