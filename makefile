CC=gcc
CFLAGS=-g -Wall -MMD -O3 -ffast-math -march=native -I/usr/X11R6/include
LDFLAGS=-g -L/usr/X11R6/lib
SOURCES=galsim.c graphics/graphics.c
OBJECTS=$(SOURCES:.c=.o)

galsim: $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@ -lm -lX11

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

-include $(SOURCES:.c=.d)

all: galsim compare

compare: compare_gal_files/compare_gal_files.c
	$(CC) -o compare_gal_files/compare_gal_files compare_gal_files/compare_gal_files.c

memcheck: galsim
	valgrind --leak-check=full ./galsim

.PHONY: clean all memcheck

clean:
	@rm -f *.o *.d galsim graphics/*.o graphics/*.d compare_gal_files/compare_gal_files