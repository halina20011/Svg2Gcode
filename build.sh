#!/bin/sh

# Create directory
mkdir -p ./Build
GCCFLAGS="-Wextra -Wall"
gcc main.c $(xml2-config --cflags --libs) $GCCFLAGS -o ./Build/svg2Gcode
