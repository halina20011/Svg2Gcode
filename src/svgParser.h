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

#ifndef SVGPARSER
#define SVGPARSER

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>       // pow, hypot
#include <float.h>      // DBL_MIN, DBL_MAX
#include <stdint.h>

#include <libxml/parser.h>

#include "./singleLinkedList.h"
#include "./vector.h"

#define BUFFER_MAX_SIZE 32

struct Geometry{
    double *points;
    size_t pointsSize;
    double bounds[4];
    double length;
};

int geometryInit(struct Geometry **geometry);
void freeGeometry(struct Geometry *geometry);
void *point2Void(double x, double y);
void calculateBounds(double *bounds, double x, double y);
void addPoint(struct Geometry *geometry, struct Vector *vector, double x, double y);

struct SVG{
    double width, height;
    uint16_t docWidth, docHeight;
    struct Geometry **geometries;
    uint16_t geometryCount;
    double bounds[4];
};

int svgInit(struct SVG **svg);
void freeSvg(struct SVG *svg);

struct Command{
    uint8_t type;
    uint8_t typeIndex;
    bool absolute;
    double *parameters;
    uint16_t parametersSize;
};

void printCommand(struct Command *command);
void freeCommand(struct Command *command);

struct Path{
    struct Command **commands;
    uint16_t commandsSize;
};

void freePath(struct Path **path);

double quadraticBezierCurve(double t, double p0, double p1, double p2);
double cubicBezierCurve(double t, double p0, double p1, double p2, double p3);
void quadraticBeziercurveto(struct Geometry *geometry, struct Vector *vector, double *x, double *y, double x1, double y1, double x2, double y2);
void curveto(struct Geometry *geometry, struct Vector *vector, double *x, double *y, double x1, double y1, double x2, double y2, double x3, double y3);
void centerParameterization(double *cx, double *cy, double *startAngle, double *deltaAngle, double x1, double y1, double x2, double y2, double rx, double ry, double phi, bool fa, bool fs);

void ellipticalArc(struct Geometry *geometry, struct Vector *vector, double *x1, double *y1, double x2, double y2, double rx, double ry, double xAxisRotation, bool largeArcFlag, bool sweepFlag);
double absolute(double currentPoint, double point, bool absolute);
void updateReflection(double *reflection, double pX3, double pY3, double pX4, double pY4);
void path2Geometry(struct Geometry ***geometries, uint16_t *geometryCount, struct Path *path);
void printPath(struct Path *path);
void parsePath(xmlNode *node, struct Path **path);
void traversElement(xmlNode *node, unsigned int depth, struct SVG *svg);

#endif
