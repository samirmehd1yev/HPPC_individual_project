# Compiler and flags
CC=/opt/homebrew/opt/llvm/bin/clang
CFLAGS=-g -Wall -MMD -O3 -ffast-math -march=native -fopenmp -I/usr/X11R6/include
LDFLAGS=-g -fopenmp -L/usr/X11R6/lib 

# Sources and objects
SOURCES=galsim.c ../graphics/graphics.c
OBJECTS=$(SOURCES:.c=.o)

# Target for galsim
galsim: $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@ -lm -lX11

# Rule for object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Include dependencies
-include $(SOURCES:.c=.d)

# Default target
all: galsim compare

# Target for compare
compare: ../compare_gal_files/compare_gal_files.c
	$(CC) -o ../compare_gal_files/compare_gal_files ../compare_gal_files/compare_gal_files.c

# Target for memory check
memcheck: galsim
	valgrind --leak-check=full ./galsim

# Phony targets
.PHONY: clean all memcheck

# Clean target
clean:
	@rm -f *.o *.d galsim ../graphics/*.o ../graphics/*.d ../compare_gal_files/compare_gal_files
