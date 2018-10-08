.PHONY: clean

all: fetch_data

fetch_data: fetch_data.c
	icc -g -O2 -Wall fetch_data.c -o fetch_data -lm

clean:
	$(RM) fetch_data *.o

