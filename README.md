# Svg2Gcode
Simple command line interface to convert svg image to gcode

## Dependencies
libxml2 </br>

## Build
```chmod +x build.sh``` </br>
```./build.sh``` </br>

## Help message
```
Arguments:
	--help            show this help message
	-f [file]         set input svg image
	-o [file]         set name of output file, if not set "output.gcode" will be used
	-w [scale]        scale image width to this size in mm
	-h [scale]        scale image height to this size in mm
	-r                keep the ratio when scalling with "-w" or "-h"
	-c                move the image to the center of given size (width and height have to be set)
	-cw [size]        specify max print width
	-ch [size]        specify max print height
	-z [offset]       set z offset when traveling in mm (default)
	-s [speed]        set speed rate for drawing in mm/min (default 1500 mm/min)
	-st [speed]       set speed rate for traveling at x, y axis in mm/min (default 5000 mm/min)
	-sz [speed]       set speed rate for traveling at z axis (up) in mm/min (default 3000 mm/min)
```

# Examples
```./Build/svg2Gcode -f Examples/loremIpsum.svg -cw 220 -ch 220 -c``` </br>
```./Build/svg2Gcode -f ./Examples/text.svg -z 10 -w 220 -h 200 -o ./Examples/example.gcode``` </br>
```./Build/svg2Gcode -f ./Examples/img1.svg``` </br>
