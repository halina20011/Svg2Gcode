// Copyright (C) 2023  halina20011
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "./singleLinkedList.h"
#include "./svgParser.h"
#include "./vector.h"

#include <libxml/parser.h>

bool centerImage = false, keepRatio = false;
double xOffset = 0.0, yOffset = 0.0, zOffset = 3.0;
uint32_t speed = 1500, tSpeed = 5000, zSpeed = 3000;
double width = -1.0, height = -1.0;
double clipWidth = -1.0, clipHeight = -1.0;

void relativeRect(double *rect, double *points){
    rect[0] = points[0];
    rect[1] = points[1];
    rect[2] = points[2] - points[0];
    rect[3] = points[3] - points[1];
}

int16_t parseInt16_t(char *strValue){
    char *endPtr;
    long val = strtol(strValue, &endPtr, 10);
    if(val < INT16_MIN || INT16_MAX < val){
        return 0;
    }

    return (int16_t)val;
}

void svgFileProperties(xmlNode *rootElement, uint16_t *docWidth, uint16_t *docHeight){
    xmlNodePtr node;
    xmlChar *widthData, *heightData;
    
    for(node = rootElement; node != NULL; node = node->next){
        if(node->type == XML_ELEMENT_NODE && xmlStrcmp(node->name, BAD_CAST"svg") == 0){
            widthData = xmlGetProp(node, (const xmlChar*)"docWidth");
            heightData = xmlGetProp(node, (const xmlChar*)"docHeight");
        }
    }
    
    // Units will be ignored
    long tempVal = 0;
    *docWidth = (widthData != NULL) ? (tempVal = parseInt16_t((char*)widthData), xmlFree(widthData), tempVal) : 0;
    *docHeight = (heightData != NULL) ? (tempVal = parseInt16_t((char*)heightData), xmlFree(heightData), tempVal) : 0;
}

void checkParameter(int argc, int i, uint8_t paramSize){
    if(argc <= i + paramSize){
        fprintf(stderr, "Wrong number of arguments after a flag\n");
        exit(1);
    }
}

void normalizeSvg(struct SVG *svg){
    for(uint16_t i = 0; i < svg->geometryCount; i++){
        struct Geometry *geometry = svg->geometries[i];
        double *bounds = geometry->bounds;
        
        if(bounds[0] < svg->bounds[0]){
            svg->bounds[0] = bounds[0];
        }
        if(svg->bounds[2] < bounds[2]){
            svg->bounds[2] = bounds[2];
        }

        if(bounds[1] < svg->bounds[1]){
            svg->bounds[1] = bounds[1];
        }
        if(svg->bounds[3] < bounds[3]){
            svg->bounds[3] = bounds[3];
        }
    }

    // printf("bounds [%g %g %g %g]\n", svg->bounds[0], svg->bounds[1], svg->bounds[2], svg->bounds[3]);
    svg->width = svg->bounds[2] - svg->bounds[0];
    svg->height = svg->bounds[3] - svg->bounds[1];

    xOffset = -svg->bounds[0];
    yOffset = -svg->bounds[1];
    svg->bounds[2] -= svg->bounds[0];
    svg->bounds[3] -= svg->bounds[1];

    svg->bounds[0] = 0.0;
    svg->bounds[1] = 0.0;

    // printf("img size %g %g\n", svg->width, svg->height);
    // printf("img relative bounds [%g %g %g %g]\n", svg->bounds[0], svg->bounds[1], svg->bounds[2], svg->bounds[3]);
    
    double xScale = 1.0, yScale = 1.0;
    if(width != -1.0 || height != -1.0){
        xScale = width / svg->width;
        yScale = height / svg->height;
        if(keepRatio){
            if(width != -1.0){
                yScale = xScale;
            }
            else{
                xScale = yScale;
            }
        }
    }
    // printf("scale [x, y] %g %g\n", xScale, yScale);

    if(centerImage){
        const double xSize = svg->width * xScale, ySize = svg->height * yScale;
        if(clipWidth < xSize || clipHeight < ySize){
            fprintf(stderr, "error: image is too big for current size cant be centered\n");
            exit(1);
        }
        
        // printf("cw %g x %g\n", clipWidth, xSize);
        // printf("cw %g x %g\n", clipHeight, ySize);
        xOffset = (clipWidth - xSize) / 2.0 + xOffset;
        yOffset = -((clipHeight - ySize) / 2.0 - yOffset);
    }
    // printf("offset [%g %g]\n", xOffset, yOffset);

    for(uint16_t g = 0; g < svg->geometryCount; g++){
        struct Geometry *geometry = svg->geometries[g];
        // printf("size %zu\n", geometry->pointsSize);
        for(size_t i = 0; i < geometry->pointsSize; i++){
            double *x = geometry->points + 2 * i + 0;
            double *y = geometry->points + 2 * i + 1;
            // printf("%zu\n", i);
            // rotate on the y axis
            *x = (*x + xOffset) * xScale;
            *y = (*y - svg->height + yOffset) * -1.0;
            *y = (*y) * yScale;
        }
    }
}

void svg2Gcode(struct SVG *svg, char *fileOutput){
    FILE *fp;
    fp = fopen(fileOutput, "w");
    if(fp == NULL){
        fprintf(stderr, "error: \"%s\" couldnt be opened\n", fileOutput);
        exit(1);
    }

    // Add file header
    fprintf(fp, "; Svg document image size: %i %i\n", svg->docWidth, svg->docHeight);
    fprintf(fp, "; Svg image size: %g %g\n", svg->width, svg->height);
    fprintf(fp, "G28; Auto home\n");
    fprintf(fp, "G90; Absolute positioning\n");
    fprintf(fp, "G21; Set units to millimeters\n");
    fprintf(fp, "M75 %s; Start print job timer\n", fileOutput);

    for(uint16_t i = 0; i < svg->geometryCount; i++){
        fprintf(fp, "; Geomtery %i\n", i);
        struct Geometry *geometry = svg->geometries[i];
        // struct Node *current = geometry->pointsHead;
        
        // for(int pointIndex = 0; current != NULL; pointIndex++){
        for(size_t pointIndex = 0; pointIndex < geometry->pointsSize; pointIndex ++){
            double x = geometry->points[2 * pointIndex + 0];
            double y = geometry->points[2 * pointIndex + 1];
            // void *data = current->data;
            // double x = **((double**) (data + 0)) + xOffset;
            // double y = **((double**) (data + sizeof(double*)));
            // printf("Points [%g %g]\n", x, y);
            // y = (y - (double)(svg->docHeight)) * -1.0 + yOffset;
            // x = (clipWidth != -1 && clipWidth <= x) ? clipWidth : x;
            // y = (clipHeight != -1 && clipHeight <= y) ? clipHeight: y;
            if(pointIndex == 0){
                fprintf(fp, "G0 Z%g F%d\n", zOffset, zSpeed);
                fprintf(fp, "G0 X%g Y%g F%d\n", x, y, tSpeed);
                fprintf(fp, "G0 Z0 F%d\n", zSpeed);
            }
            else{
                fprintf(fp, "G1 X%g Y%g F%d\n", x, y, speed);
            }
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "M76; Pause print job timer\n");
    fprintf(fp, "; END\n");
    fprintf(fp, "G0 Z%g\n", zOffset);
    fprintf(fp, "G0 X0 Y%g\n", ((clipHeight != -1.0) ? clipHeight : 0.0));
    
    fclose(fp);
}

void printHelp(){
    printf("Usage: svg2Gcode [ARGUMENT]\n");
    printf("Arguments:\n");
    printf("\t--help            show this help message\n");
    printf("\t-f [file]         set input svg image\n");
    printf("\t-o [file]         set name of output file, if not set \"output.gcode\" will be used\n");

    printf("\t-w [scale]        scale image width to this size in mm\n");
    printf("\t-h [scale]        scale image height to this size in mm\n");
    printf("\t-r                keep the ratio when scalling with \"-w\" or \"-h\"\n");

    printf("\t-c                move the image to the center of given size (width and height have to be set)\n");

    printf("\t-cw [size]        specify max print width\n");
    printf("\t-ch [size]        specify max print height\n");

    printf("\t-z [offset]       set z offset when traveling in mm (default)\n");
    printf("\t-s [speed]        set speed rate for drawing in mm/min (default %d mm/min)\n", speed);
    printf("\t-st [speed]       set speed rate for traveling at x, y axis in mm/min (default %d mm/min)\n", tSpeed);
    printf("\t-sz [speed]       set speed rate for traveling at z axis (up) in mm/min (default %d mm/min)\n", zSpeed);
}

int main(int argc, char **argv){
    // uint8_t usePadding = 0;
    // double padding[4] = {-1.0};
        
    char *fileName = NULL;
    char *fileOutput = NULL;
    if(argc < 2){
        fprintf(stderr, "Wrong number of arguments\n");
        return 1;
    }

    for(int i = 1; i < argc; i++){
        if(argv[i][0] != '-'){
            fprintf(stderr, "Invalid flag format\n");
            return 1;
        }
        if(strcmp(argv[i], "--help") == 0){
            printHelp();
            return 0;
        }
        if(strcmp(argv[i], "-f") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            fileName = argv[i];
        }
        else if(strcmp(argv[i], "-o") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            fileOutput = argv[i];
        }
        else if(strcmp(argv[i], "-w") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            width = strtod(argv[i], NULL);
        }
        else if(strcmp(argv[i], "-h") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            height = strtod(argv[i], NULL);
        }
        else if(strcmp(argv[i], "-r") == 0){
            keepRatio = true;
        }
        else if(strcmp(argv[i], "-cw") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            clipWidth = strtod(argv[i], NULL);
        }
        else if(strcmp(argv[i], "-ch") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            clipHeight = strtod(argv[i], NULL);
        }
        else if(strcmp(argv[i], "-z") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            zOffset = strtod(argv[i], NULL);
        }
        else if(strcmp(argv[i], "-c") == 0){
            centerImage = true;
        }
        else if(strcmp(argv[i], "-s") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            speed = strtol(argv[i], NULL, 10);
        }
        else if(strcmp(argv[i], "-st") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            tSpeed = strtol(argv[i], NULL, 10);
        }
        else if(strcmp(argv[i], "-sz") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            zSpeed = strtol(argv[i], NULL, 10);
        }
        else{
            fprintf(stderr, "invalid option: \"%s\"\n", argv[i]);
            fprintf(stderr, "Try: '%s --help' for more info\n", argv[0]);
            return 1;
        }
        // else if(strncmp(argv[i], "-p", 2) == 0){
        //     checkParameter(argc, i, 4);
        //     usePadding = (strncmp(argv[i], "-pa", 2) == 0) ? 2 : 1;
        //     printBounds[0] = strtod(argv[i + 1], NULL);
        //     printBounds[1] = strtod(argv[i + 2], NULL);
        //     printBounds[2] = strtod(argv[i + 3], NULL);
        //     printBounds[3] = strtod(argv[i + 4], NULL);
        //     i += 4;
        // }
    }

    if(fileName == NULL){
        fprintf(stderr, "Missing input file\n");
        return 1;
    }

    if(fileOutput == NULL){
        fileOutput = "output.gcode";
    }

    xmlDoc *document = NULL;
    xmlNode *rootElement = NULL;
    
    document = xmlReadFile(fileName, NULL, 0);
    if(document == NULL){
        fprintf(stderr, "error: couldnt parse file %s\n", fileName);
        return 1;
    }

    rootElement = xmlDocGetRootElement(document);

    struct SVG *svg = NULL;
    if(svgInit(&svg)){
        fprintf(stderr, "error: new svg couldnt be created\n");
        return 1;
    }

    if(centerImage){
        if(clipWidth == -1.0 || clipHeight == -1.0){
            fprintf(stderr, "error: clipWidth and clipHeight have to be set when using center image\n");
            return 1;
        }
    }

    svgFileProperties(rootElement, &svg->docWidth, &svg->docHeight);
    
    // result will be and svg struct that has an array of geometries =>
    // geometry is list of points that will have start in left up corner
    // reverse the geometry on y axis
    traversElement(rootElement, 0, svg);
    normalizeSvg(svg);
    // printf("[%g %g %g %g]\n", svg->bounds[0], svg->bounds[1], svg->bounds[2], svg->bounds[3]);
    // transformSvg(svg, width, height, usePadding, padding);

    svg2Gcode(svg, fileOutput);

    freeSvg(svg);
    xmlFreeDoc(document);
    xmlCleanupParser();

    return 0;
}
