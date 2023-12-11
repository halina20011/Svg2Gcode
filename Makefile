CC = gcc
CFLAGS = -Wall -Wextra -Wshadow -Isrc
LIBS = $(shell xml2-config --cflags --libs) -lm

FILES = $(wildcard ./src/*.c)
HEADERS = $(patsubst ./src/%.c, ./build/%.h, $(FILES))
OBJECTS = $(patsubst ./src/%.c, ./build/%.o, $(FILES))

debug:
	@echo "files $(FILES)"
	@echo "objects $(OBJECTS)"

buildDir:
	mkdir -p build

build/%.o: src/%.c buildDir
	echo "building $@"
	$(CC) $(CFLAGS) $(LIBS) -c -o $@ $<

main: $(FILES) $(OBJECTS)
	echo $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) $(LIBS) -o ./build/svg2Gcode

build: main
	echo $(FILES)
	echo $(OBJECTS)
