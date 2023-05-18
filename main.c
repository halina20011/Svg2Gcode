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

#include <libxml/parser.h>

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

void svgFileProperties(xmlNode *rootElement, uint16_t *width, uint16_t *height){
    xmlNodePtr node;
    xmlChar *widthData, *heightData;
    
    for(node = rootElement; node != NULL; node = node->next){
        if(node->type == XML_ELEMENT_NODE && xmlStrcmp(node->name, BAD_CAST"svg") == 0){
            widthData = xmlGetProp(node, (const xmlChar*)"width");
            heightData = xmlGetProp(node, (const xmlChar*)"height");
        }
    }
    
    // Units will be ignored
    long tempVal = 0;
    *width = (widthData != NULL) ? (tempVal = parseInt16_t((char*)widthData), xmlFree(widthData), tempVal) : 0;
    *height = (heightData != NULL) ? (tempVal = parseInt16_t((char*)heightData), xmlFree(heightData), tempVal) : 0;
}

void checkParameter(int argc, int i, uint8_t paramSize){
    if(argc <= i + paramSize){
        fprintf(stderr, "Wrong number of arguments after a flag\n");
        exit(1);
    }
}

void svg2Gcode(struct SVG *svg, char *fileOutput, double width, double height, double zOffset){
    FILE *fp;
    fp = fopen(fileOutput, "w");
    if(fp == NULL){
        exit(1);
    }

    // Add file header
    fprintf(fp, "; Svg image size: %i %i\n", svg->width, svg->height);
    fprintf(fp, "G28; Auto home\n");
    fprintf(fp, "G90; Absolute positioning\n");
    fprintf(fp, "G21; Set units to millimeters\n");
    fprintf(fp, "M75 %s; Start print job timer\n", fileOutput);

    for(uint16_t i = 0; i < svg->geometryCount; i++){
        fprintf(fp, "; Geomtery %i\n", i);
        struct Geometry *geometry = svg->geometries[i];
        Node *current = geometry->pointsHead;

        for(int pointIndex = 0; current->next != NULL; pointIndex++){
            current = current->next;
            void *data = current->data;
            double x = **((double**) (data + 0));
            double y = **((double**) (data + sizeof(double*)));
            // printf("Points [%g %g]\n", x, y);
            y = (y - svg->height) * -1.0;
            x = (width != -1 && width <= x) ? width : x;
            y = (height != -1 && height <= y) ? height: y;
            if(pointIndex == 0){
                fprintf(fp, "G0 Z%g\n", zOffset);
                fprintf(fp, "G0 X%g Y%g\n", x, y);
                fprintf(fp, "G0 Z0\n");
            }
            fprintf(fp, "G1 X%g Y%g\n", x, y);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "M76; Pause print job timer\n");
    fprintf(fp, "; END\n");
    fprintf(fp, "G0 Z%g\n", zOffset);
    fprintf(fp, "G0 X0 Y%g\n", ((height != -1) ? height : 0.0));
    
    fclose(fp);
}

void printHelp(){
    printf("Usage: svg2Gcode [ARGUMENT]\n");
    printf("Arguments:\n");
    printf("\t--help            show this help message\n");
    printf("\t-f [file]         set input svg image\n");
    printf("\t-o [file]         set name of output file, if not set \"output.gcode\" will be used\n");
    printf("\t-w [size]         specify max print width\n");
    printf("\t-h [size]         specify max print height\n");
    printf("\t-z [offset]       set z offset when traveling\n");
}

int main(int argc, char **argv){
    // printf("I live :3\n");
    
    // uint8_t usePadding = 0;
    // double padding[4] = {-1.0};
    double width = -1.0;
    double height = -1.0;
    double zOffset = 5.0;

    char *fileName = NULL;
    char *fileOutput = NULL;
    char *end;
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
            width = strtod(argv[i], &end);
        }
        else if(strcmp(argv[i], "-h") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            height = strtod(argv[i], &end);
        }
        else if(strcmp(argv[i], "-z") == 0){
            checkParameter(argc, i, 1);
            i += 1;
            zOffset = strtod(argv[i], &end);
        }
        else{
            fprintf(stderr, "invalid option: \"%s\"\n", argv[i]);
            fprintf(stderr, "Try: '%s --help' for more info\n", argv[0]);
            return 1;
        }
        // else if(strncmp(argv[i], "-p", 2) == 0){
        //     checkParameter(argc, i, 4);
        //     usePadding = (strncmp(argv[i], "-pa", 2) == 0) ? 2 : 1;
        //     printBounds[0] = strtod(argv[i + 1], &end);
        //     printBounds[1] = strtod(argv[i + 2], &end);
        //     printBounds[2] = strtod(argv[i + 3], &end);
        //     printBounds[3] = strtod(argv[i + 4], &end);
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
        fprintf(stderr, "error: couldnt create new svg\n");
        return 1;
    }

    svgFileProperties(rootElement, &svg->width, &svg->height);

    traversElement(rootElement, 0, svg);
    updateSvgBound(svg);
    // transformSvg(svg, width, height, usePadding, padding);

    // double bounds[4] = {0};
    // relativeRect(bounds, svg->bounds);
    // printf("[%g %g %g %g]\n", bounds[0], bounds[1], bounds[2], bounds[3]);

    svg2Gcode(svg, fileOutput, width, height, zOffset);

    freeSvg(svg);
    xmlFreeDoc(document);
    xmlCleanupParser();

    return 0;
}
