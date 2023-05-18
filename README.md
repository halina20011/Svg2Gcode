# Svg2Gcode
Simple command line interface to convert svg image to gcode

## Dependencies
libxml2 </br>

## Build
```chmod +x build.sh``` </br>
```./build.sh``` </br>

## Help message
```
Usage: svg2Gcode [ARGUMENT]
Arguments:
	--help            show this help message
	-f [file]         set input svg image
	-o [file]         set name of output file, if not set "output.gcode" will be used
	-w [size]         specify max print width
	-h [size]         specify max print height
	-z [offset]       set z offset when traveling
```

# Examples
```./Build/svg2Gcode -f ./Examples/text.svg -z 10 -w 220 -h 200 -o ./Examples/example.gcode``` </br>
```./Build/svg2Gcode -f ./Examples/img1.svg``` </br>
