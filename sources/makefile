CC=g++

all: genesplicer

genesplicer: genesplicer.o sites.o
	$(CC) -o genesplicer genesplicer.o sites.o -lm

genesplicer.o: genesplicer.cpp
	$(CC) -c genesplicer.cpp

sites.o: sites.c
	$(CC) -c sites.c
