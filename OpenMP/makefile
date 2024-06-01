CC=gcc
CFLAGS=-g -Wall -MMD -O3 -ffast-math -march=native -fopenmp
LDFLAGS=-g -fopenmp

SOURCES=galsim.c ../graphics/graphics.c
OBJECTS=$(SOURCES:.c=.o)

galsim: $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@ -lm -lX11

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

-include $(SOURCES:.c=.d)

all: galsim compare

compare: ../compare_gal_files/compare_gal_files.c
	$(CC) -o ../compare_gal_files/compare_gal_files ../compare_gal_files/compare_gal_files.c -lm

memcheck: galsim
	valgrind --leak-check=full ./galsim

.PHONY: clean all memcheck

clean:
	@rm -f *.o *.d galsim ../graphics/*.o ../graphics/*.d ../compare_gal_files/compare_gal_files
