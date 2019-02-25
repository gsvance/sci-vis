.PHONY: clean

all: get_sdf_data

get_sdf_data: get_sdf_data.c
	icc -g -O2 -Wall get_sdf_data.c -o get_sdf_data -lm

clean:
	$(RM) get_sdf_data *.o

