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

#include "svgParser.h"
#include "mathFunc.h"

#define COMMANDS_LENGHT 10
uint8_t commandsTypes[COMMANDS_LENGHT]  = {'m', 'z', 'l', 'h', 'v', 'c', 's', 'q', 't', 'a'};
uint8_t commandsSizes[COMMANDS_LENGHT]  = {  2,   0,   2,   1,   1,   6,   4,   4,   2,   7};

int geometryInit(struct Geometry **geometry){
    *geometry = malloc(sizeof(struct Geometry));
    if(*geometry == NULL){
        return 1;
    }
    
    (*geometry)->pointsSize = 0;
    (*geometry)->bounds[0] = DBL_MAX;
    (*geometry)->bounds[1] = DBL_MAX;
    (*geometry)->bounds[2] = DBL_MIN;
    (*geometry)->bounds[3] = DBL_MIN;
    (*geometry)->length =  0;

    return 0;
}

void freeGeometry(struct Geometry *geometry){
    free(geometry->points);
}

void *point2Void(double x, double y){
    double *newX = malloc(sizeof(double));
    double *newY = malloc(sizeof(double));
    if(newX == NULL || newY == NULL){
        return NULL;
    }
    *newX = x;
    *newY = y;

    void *copy = malloc(sizeof(double*) * 2);
    *((double**) (copy + 0)) = newX;
    *((double**) (copy + sizeof(double))) = newY;

    return copy;
}

//     3
//    ---
// 0 |   | 2
//    ---
//     1
void calculateBounds(double *bounds, double x, double y){
    if(x < bounds[0]){
        bounds[0] = x;
    }
    if(bounds[2] < x){
        bounds[2] = x;
    }

    if(y < bounds[1]){
        bounds[1] = y;
    }
    if(bounds[3] < y){
        bounds[3] = y;
    }
}

void addPoint(struct Geometry *geometry, struct Vector *vector, double x, double y){
    calculateBounds(geometry->bounds, x, y);
    double d = hypot(x, y);

    geometry->length += d;
    
    vectorPush(vector, x);
    vectorPush(vector, y);
    // void *voidPoint = point2Void(x, y);
    // printf("adding [%g %g] point\n", x, y);
    // push(&geometry->pointsHead, &geometry->pointsTail, voidPoint);
}

int svgInit(struct SVG **svg){
    *svg = malloc(sizeof(struct SVG));
    if(*svg == NULL){
        return 1;
    }
    (*svg)->width = 0;
    (*svg)->height = 0;
    (*svg)->geometryCount = 0;
    (*svg)->geometries = malloc(0);

    (*svg)->bounds[0] = DBL_MAX;
    (*svg)->bounds[1] = DBL_MAX;
    (*svg)->bounds[2] = DBL_MIN;
    (*svg)->bounds[3] = DBL_MIN;

    return 0;
}

void freeSvg(struct SVG *svg){
    for(uint16_t i = 0; i < svg->geometryCount; i++){
        freeGeometry(svg->geometries[i]);
    }
    free(svg);
}

// void removePaddingSvg(struct SVG *svg){
//     double offsetX = svg->bounds[0];
//     double offsetY = svg->bounds[1];
//     double sizeX = svg->bounds[2] - svg->bounds[0];
//     double sizeY = svg->bounds[3] - svg->bounds[1];
// }

// Scale the image
// void transformSvg(struct SVG *svg, double width, double height, bool usePadding, double *padding){
//     // if its set -p (usePadding 1) then set new padding by re calculating the new the offset of x and y
//     // if width and height are not set then just flip the y axis so 0,0 will be in left bottom corner
//     // if only width is set then scale the image to the width and keep ratio
//
//     double offsetX = svg->bounds[0];
//     // double scale = 1.0 / ((sizeY < sizeX) ? sizeX : sizeY);
//     double scale = 1.0 / (double)((svg->height < svg->width) ? svg->width : svg->height);
//     printf("S: %g\n", scale);
//     for(uint16_t i = 0; i < svg->geometryCount; i++){
//         struct Geometry *geometry = svg->geometries[i];
//         struct Node *current = geometry->pointsHead;
//         while(current->next != NULL){
//             current = current->next;
//             void *data = current->data;
//             double *x = *((double**) (data + 0));
//             double *y = *((double**) (data + sizeof(double*)));
//         }
//     }
// }

void printCommand(struct Command *command){
    printf("%c ", command->type);
    for(int i = 0; i < command->parametersSize; i++){
        printf("%g ", command->parameters[i]);
    }

    printf("\n");
}

void freeCommand(struct Command *command){
    free(command->parameters);
    free(command);
}

void freePath(struct Path **path){
    for(uint16_t i = 0; i < (*path)->commandsSize; i++){
        freeCommand((*path)->commands[i]);
    }

    free((*path)->commands);
    free(*path);
}

double quadraticBezierCurve(double t, double p0, double p1, double p2){
    double point0 = pow((1.0 - t), 2.0) * p0;
    double point1 = 2.0 * (1.0 - t) * t * p1;
    double point2 = pow(t, 2.0) * p2;

    return point0 + point1 + point2;
}

double cubicBezierCurve(double t, double p0, double p1, double p2, double p3){
    double point0 = pow((1.0 - t), 3.0) * p0;
    double point1 = 3.0 * pow(1.0 - t, 2.0) * t * p1;
    double point2 = 3.0 * (1.0 - t) * pow(t, 2.0) * p2;
    double point3 = pow(t, 3.0) * p3;

    return point0 + point1 + point2 + point3;
}

// TODO calculate curveto resolution
double resolution = 0.1;

void quadraticBeziercurveto(struct Geometry *geometry, struct Vector *vector, double *x, double *y, double x1, double y1, double x2, double y2){
    for(double t = 0.0; t <= 1.0 + resolution; t = fmin(t + resolution, 1.0)){
        double px = quadraticBezierCurve(t, *x, x1, x2);
        double py = quadraticBezierCurve(t, *y, y1, y2);
        addPoint(geometry, vector, px, py);
        if(1.0 <= t){
            *x = px;
            *y = py;
            break;
        }
    }
}

void curveto(struct Geometry *geometry, struct Vector *vector, double *x, double *y, double x1, double y1, double x2, double y2, double x3, double y3){
    for(double t = 0.0; t <= 1.0 + resolution; t = fmin(t + resolution, 1.0)){
        double px = cubicBezierCurve(t, *x, x1, x2, x3);
        double py = cubicBezierCurve(t, *y, y1, y2, y3);
        addPoint(geometry, vector, px, py);
        if(1.0 <= t){
            *x = px;
            *y = py;
            break;
        }
    }
}

// https://www.w3.org/TR/SVG11/implnote.html#ArcConversionEndpointToCenter
void centerParameterization(double *cx, double *cy, double *startAngle, double *deltaAngle, double x1, double y1, double x2, double y2, double rx, double ry, double phi, bool fa, bool fs){
    // printf("%g, %g, %g, %g, %g, %g, %g, %i, %i\n", x1, y1, x2, y2, rx, ry, phi, (int)fa, (int)fs);
    double phiSin = sin(phi);
    double phiCos = cos(phi);

    double xHalfDiff = (x1 - x2) / 2.0;
    double yHalfDiff = (y1 - y2) / 2.0;
    double xHalfSum  = (x1 + x2) / 2.0;
    double yHalfSum  = (y1 + y2) / 2.0;

    // Compute x1', y1'
    double x1Prime = phiCos * xHalfDiff + phiSin * yHalfDiff;
    double y1Prime = -phiSin * xHalfDiff + phiCos * yHalfDiff;
    // printf("%g %g\n", x1Prime, y1Prime);

    double x1PrimeSquare = pow(x1Prime, 2);
    double y1PrimeSquare = pow(y1Prime, 2);

    // Radii out of range correction
    rx = fabs(rx);
    ry = fabs(ry);

    double lambda = (x1PrimeSquare / (rx*rx)) + (y1PrimeSquare / (ry*ry));
    if(1.0 < lambda){
        double lambdaSqrt = sqrt(lambda);
        rx = lambdaSqrt * rx;
        ry = lambdaSqrt * ry;
    }

    double rxSquare = pow(rx, 2);
    double rySquare = pow(ry, 2);

    // Compute cx', cy'
    double line1 = (rxSquare * rySquare) - (rxSquare * y1PrimeSquare) - (rySquare * x1PrimeSquare);
    double line2 = (rxSquare * y1PrimeSquare) + (rySquare * x1PrimeSquare);
    double scalar = sqrt(fabs(line1/line2));
    if(fa == fs){
        scalar *= -1.0;
    }

    // printf("L %g %g Scalar: %g\n", line1, line2, scalar);

    double cxPrime = scalar * ((rx * y1Prime) / ry);
    double cyPrime = scalar * (-(ry * x1Prime) / rx);
    // printf("%g %g\n", cxPrime, cyPrime);

    *cx = (phiCos * cxPrime) + (-phiSin * cyPrime) + xHalfSum;
    *cy = (phiSin * cxPrime) + (phiCos * cyPrime) + yHalfSum;

    // Calculate start and end angle
    double xVector =  (x1Prime - cxPrime) / rx;
    double yVector =  (y1Prime - cyPrime) / ry;
    double xVector2 = (x1Prime + cxPrime) / rx;
    double yVector2 = (y1Prime + cyPrime) / ry;
    // printf("%g %g %g %g\n", xVector, yVector, xVector2, yVector2);
    
    *startAngle = rad(angleCalc(1.0, 0.0, xVector, yVector));
    *deltaAngle = rad(angleCalc(xVector, yVector, -xVector2, -yVector2));
    normalizeAngle360(deltaAngle);
    if(fs == false){
        *deltaAngle -= 360.0;
    }

    // printf("Center of ellipse is [%g %g] %g %g\n", *cx, *cy, *startAngle, *deltaAngle);
}

void ellipticalArc(struct Geometry *geometry, struct Vector *vector, double *x1, double *y1, double x2, double y2, double rx, double ry, double xAxisRotation, bool largeArcFlag, bool sweepFlag){
    double cx, cy, startAngle, deltaAngle;
    // Calculate center point
    centerParameterization(&cx, &cy, &startAngle, &deltaAngle, *x1, *y1, x2, y2, rx, ry, xAxisRotation, largeArcFlag, sweepFlag);
    
    double endAngle = startAngle + deltaAngle;
    normalizeAngle360(&endAngle);
    // printf("Angles: %g %g %g\n", startAngle, deltaAngle, endAngle);

    double theta = xAxisRotation;
    // TODO: precalculate resolution of elliptical arc
    double re = 2.0;
    bool clockwise = sweepFlag;
    re *= (clockwise) ? 1.0 : -1.0;
    double angle = 0.0;
    double (*compareFunc)(double, double) = (clockwise) ? &fmin : &fmax;
    while(((clockwise) ? (angle < deltaAngle + re) : (deltaAngle - re < angle))){
        double alpha = ((startAngle + angle) * M_PI) / 180;
        double xr = rx * cos(alpha);
        double yr = ry * sin(alpha);

        double x = cos(theta) * xr + -sin(theta) * yr + cx;
        double y = sin(theta) * xr + cos(theta) * yr + cy;

        addPoint(geometry, vector, x, y);
        
        if(deltaAngle == angle){
            break;
        }

        angle = compareFunc(angle + re, deltaAngle);
    }

    *x1 = x2;
    *y1 = y2;
}

double absolute(double currentPoint, double point, bool absolute){
    return (absolute == true) ? point : currentPoint + point;
}

void updateReflection(double *reflection, double pX3, double pY3, double pX4, double pY4){
    reflection[0] = pX4 - pX3;
    reflection[1] = pY4 - pY3;
}

// loop throw all commands and convert them to geometry
// the geometry is closed when there is a z command
void path2Geometry(struct Geometry ***geometries, uint16_t *geometryCount, struct Path *path){
    double x = 0, y = 0;
    double commandStartX = 0, commandStartY = 0;

    char lastCommandType;
    double curvetoReflection[2] = {0};
    double bezierCurveReflection[2] = {0};

    struct Geometry *geometry;
    struct Vector *vector = NULL;
    geometryInit(&geometry);
    newVector(&vector);
    
    uint16_t size = 1;

    struct Command *command = NULL;
    struct Node *geometriesHead = NULL, *geometriesTail = NULL;

    for(uint16_t commandIndex = 0; commandIndex < path->commandsSize; commandIndex++){
        command = path->commands[commandIndex];
        // printCommand(command);
        uint16_t parametersSize = command->parametersSize;
        uint8_t commandSize = commandsSizes[command->typeIndex];
        // printf("typeIndex: %u\n", command->typeIndex);
        // printf("Command[%u]: %c\n", parametersSize, (char)command->type);

        if(commandSize != 0 && parametersSize % commandSize != 0){
            fprintf(stderr, "%c command has wrong number of parameters\n", (char)command->type);
            exit(1);
        }

        if(command->type == 'm'){
            // printf("M command\n");
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                x = command->parameters[i + 0] + ((command->absolute) ? 0.0 : x);
                y = command->parameters[i + 1] + ((command->absolute) ? 0.0 : y);
                if(i == 0){
                    commandStartX = x;
                    commandStartY = y;
                }
                addPoint(geometry, vector, x, y);
                // printf("Moved to [%g %g]\n", x, y);
            }
        }
        else if(command->type == 'z'){
            // printf("closing path on [%g %g]\n", commandStartX, commandStartY);
            addPoint(geometry, vector, commandStartX, commandStartY);
                
            if(commandIndex != path->commandsSize - 1){
                // printf("new geometry alloc\n");

                // void *voidPoint = point2Void(commandStartX, commandStartY);
                // push(&geometry->pointsHead, &geometry->pointsTail, voidPoint);
                geometry->points = vectorGetDataAndFree(vector, &geometry->pointsSize);
                push(&geometriesHead, &geometriesTail, (void*)geometry);

                geometryInit(&geometry);
                newVector(&vector);

                size += 1;
                // *geometries = realloc(*geometries, sizeof(struct Geometry*) * (*geometryCount));

                if(commandIndex + 1 < path->commandsSize && path->commands[commandIndex + 1]->type == 'm'){ 
                    x = commandStartX;
                    y = commandStartY;
                }
                else{
                    addPoint(geometry, vector, commandStartX, commandStartY);
                    x = 0.0;
                    y = 0.0;
                }
            }
        }
        else if(command->type == 'l'){
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                // printf("abs: %i\n", command->absolute);
                x = command->parameters[i + 0] + ((command->absolute) ? 0.0 : x);
                y = command->parameters[i + 1] + ((command->absolute) ? 0.0 : y);

                addPoint(geometry, vector, x, y);
                // printf("lineto [%g %g]\n", x, y);
            }
        }
        else if(command->type == 'h'){
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                x = command->parameters[i] + ((command->absolute) ? 0.0 : x);

                addPoint(geometry, vector, x, y);
                // printf("horizontal line to [%g]\n", x);
            }
        }
        else if(command->type == 'v'){
            // printf("vertical lineto\n");
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                y = command->parameters[i] + ((command->absolute) ? 0.0 : y);

                addPoint(geometry, vector, x, y);
                // printf("vertical line to [%g]\n", y);
            }
        }
        else if(command->type == 'c'){
            // printf("curveto\n");
            double *points = command->parameters;
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                double pX2 = absolute(x, points[i + 0], command->absolute);
                double pY2 = absolute(y, points[i + 1], command->absolute);
                double pX3 = absolute(x, points[i + 2], command->absolute);
                double pY3 = absolute(y, points[i + 3], command->absolute);
                double pX4 = absolute(x, points[i + 4], command->absolute);
                double pY4 = absolute(y, points[i + 5], command->absolute);
                updateReflection(curvetoReflection, pX3, pY3, pX4, pY4);
                curveto(geometry, vector, &x, &y, pX2, pY2, pX3, pY3, pX4, pY4);
                // printf("Curveto: %g %g\n", x, y);
            }
        }
        else if(command->type == 's'){
            double *points = command->parameters;
            double pX2, pY2, pX3, pY3, pX4, pY4;
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                // If command before was curveto or shorthand/smooth curveto then only use the calculated curvetoReflection
                // else the first control point will be coincident with current point (x, y)
                if((i == 0 && (lastCommandType == 'c' || lastCommandType == 's')) || 0 < i){
                    pX2 = curvetoReflection[0] + x;
                    pY2 = curvetoReflection[1] + y;
                }
                else{
                    pX2 = x;
                    pY2 = y;
                }
                pX3 = absolute(x, points[i + 0], command->absolute);
                pY3 = absolute(y, points[i + 1], command->absolute);
                pX4 = absolute(x, points[i + 2], command->absolute);
                pY4 = absolute(y, points[i + 3], command->absolute);
                updateReflection(curvetoReflection, pX3, pY3, pX4, pY4);
                curveto(geometry, vector, &x, &y, pX2, pY2, pX3, pY3, pX4, pY4);
            }
        }
        else if(command->type == 'q'){
            double *points = command->parameters;
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                double pX2 = absolute(x, points[i + 0], command->absolute);
                double pY2 = absolute(y, points[i + 1], command->absolute);
                double pX3 = absolute(x, points[i + 2], command->absolute);
                double pY3 = absolute(y, points[i + 3], command->absolute);
                updateReflection(bezierCurveReflection, pX2, pY2, pX3, pY3);
                quadraticBeziercurveto(geometry, vector, &x, &y, pX2, pY2, pX3, pY3);
            }
        }
        else if(command->type == 't'){
            double *points = command->parameters;
            double pX2, pY2, pX3, pY3;
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                // If command before was quadratic bezier curve or shorthand/quadratic bezier curve then only use the calculated bezierCurveReflection
                // else the first control point will be coincident with current point (x, y)
                if((i == 0 && (lastCommandType == 'q' || lastCommandType == 't')) || 0 < i){
                    pX2 = bezierCurveReflection[0] + x;
                    pY2 = bezierCurveReflection[1] + y;
                }
                else{
                    pX2 = x;
                    pY2 = y;
                }
                pX3 = absolute(x, points[i + 0], command->absolute);
                pY3 = absolute(y, points[i + 1], command->absolute);
                updateReflection(bezierCurveReflection, pX2, pY2, pX3, pY3);
                quadraticBeziercurveto(geometry, vector, &x, &y, pX2, pY2, pX3, pY3);
            }
        }
        else if(command->type == 'a'){
            double *points = command->parameters;
            for(uint16_t i = 0; i < parametersSize; i += commandSize){
                double rx = points[i + 0];
                double ry = points[i + 1];
                double xAxisRotation = points[i + 2];
                bool largeArcFlag = (points[i + 3] < 1) ? false : true;
                bool sweepFlag    = (points[i + 4] < 1) ?  false : true;
                double x2 = absolute(x, points[i + 5], command->absolute);
                double y2 = absolute(y, points[i + 6], command->absolute);
                // printf("Elliptical arc %g, %g, %g, %g, %g, %g, %g, %i, %i\n", x, y, x2, y2, rx, ry, xAxisRotation, (int)largeArcFlag, (int)sweepFlag);
                ellipticalArc(geometry, vector, &x, &y, x2, y2, rx, ry, xAxisRotation, largeArcFlag, sweepFlag);
            }
        }
        else{
            printf("Command \"%c\" not found\n", command->type);
        }

        // printf("[%g %g] ", y, x);
        lastCommandType = command->type;
    }

    geometry->points = vectorGetDataAndFree(vector, &geometry->pointsSize);
    push(&geometriesHead, &geometriesTail, (void*)geometry);
    
    *geometries = malloc(sizeof(struct Geometry*) * size);
    struct Node *curr = geometriesHead, *temp = NULL;
    for(int i = 0; curr; i++){
        temp = curr;
        curr = curr->next;
        (*geometries)[i] = (struct Geometry*)(temp->data);
        // free(temp);
    }

    *geometryCount = size;
}

void printPath(struct Path *path){
    for(int i = 0; i < path->commandsSize; i++){
        printCommand(path->commands[i]);
    }
}

// count number of commands in path data
// split them into seperate strings
void parsePath(xmlNode *node, struct Path **path){
    char *pathData = (char*)xmlGetProp(node, (const xmlChar*)"d");
    
    struct Node *head = NULL, *tail = NULL;

    char *buffer = malloc(sizeof(char) * (BUFFER_MAX_SIZE + 1));
    uint8_t bufferSize = 0;
    uint32_t commandsCount = 0;

    char *ptr = pathData;
    while(true){
        if(*ptr != '\0'){
            bool upperCaseChar = (65 <= (int)*ptr && (int)*ptr <= 90);
            bool lowerCaseChar = (97 <= (int)*ptr && (int)*ptr <= 122);
            bool isChar = (upperCaseChar || lowerCaseChar);

            if((isChar && *ptr != 'e') || *ptr == ' ' || *ptr == ','){
                if(0 < bufferSize){
                    buffer[bufferSize] = '\0';
                    push(&head, &tail, buffer);
                    buffer = NULL;
                    bufferSize = 0;
                    buffer = malloc(sizeof(char) * (BUFFER_MAX_SIZE + 1));
                }

                if(isChar){
                    char *command = malloc(sizeof(char) * (1 + 1));
                    command[0] = *ptr;
                    command[1] = '\0';
                    push(&head, &tail, command);
                    // printf("command added %s\n", command);
                    commandsCount++;
                }
            }
            else{
                if(bufferSize < BUFFER_MAX_SIZE){
                    buffer[bufferSize++] = *ptr;
                }
            }

            ptr += 1;
        }
        else{
            if(0 < bufferSize){
                buffer[bufferSize] = '\0';
                push(&head, &tail, buffer);
            }
            break;
        }
    }

    // printSingleLinkedList(head);

    // printf("Number of commands: %u\n", commandsCount);
    struct Command **commands = malloc(sizeof(struct Command*) * commandsCount);

    // Index thats maps the command to other information about it
    uint8_t commandTypeIndex = 0;
    uint16_t parametersIndex = 0;
    uint16_t parametersSize = 0;

    char *deletedNodesData = NULL;
    for(int i = -1; (deletedNodesData = deleteFront(&head, &tail)) != NULL; ){
        uint8_t charV = (int8_t)deletedNodesData[0];
        bool upperCaseChar = (65 <= charV && charV <= 90);
        bool lowerCaseChar = (97 <= charV && charV <= 122);
        bool isChar = (upperCaseChar || lowerCaseChar);
        charV = (upperCaseChar) ? charV + 32 : charV;

        // If its a char => new command starts
        // first get type of the command
        // make new command and set a information aboout it (is it absolute)
        // allocate first minimul size of memory that the command needs
        if(isChar){
            // printf("New command[%i] %c\n", i, (char)charV);
            for(uint8_t tI = 0; tI < COMMANDS_LENGHT; tI++){
                if(charV == commandsTypes[tI]){
                    // printf("Char %c\n", charV);
                    commandTypeIndex = tI;
                    break;
                }
            }

            commands[++i] = malloc(sizeof(struct Command));
            
            parametersSize = commandsSizes[commandTypeIndex];
            parametersIndex = 0;
            commands[i]->parameters = malloc(sizeof(double) * parametersSize);
            commands[i]->parametersSize = parametersSize;
            commands[i]->type = charV;
            commands[i]->typeIndex = commandTypeIndex;
            commands[i]->absolute = upperCaseChar;
        }
        else{
            // There is a posibility that a command will have more repeating parametersSize (x y)+
            // so if its bigger then add another layer of parameters
            if(parametersSize <= parametersIndex){
                parametersSize += commandsSizes[commandTypeIndex];
                commands[i]->parameters = realloc(commands[i]->parameters, sizeof(double) * parametersSize);
                commands[i]->parametersSize = parametersSize;
                // printf("Realocating parameters to %u\n", parametersSize);
            }

            // Just parse the char to double
            char *strEnd;
            double parameter = strtod(deletedNodesData, &strEnd);
            // printf("New parsed parameter with value %g and index %u\n", parameter, parametersIndex);
            commands[i]->parameters[parametersIndex++] = parameter;
        }

        // printf("%s\n", deletedNodesData);
        free(deletedNodesData);
    }
    
    freeSingleLinkedList(head);

    *path = malloc(sizeof(struct Path));
    (*path)->commands = commands;
    (*path)->commandsSize = commandsCount;

    // printPath(*path);

    if(pathData != NULL){
        free(pathData);
    }
}

void traversElement(xmlNode *node, unsigned int depth, struct SVG *svg){
    xmlNode *currentNode = NULL;

    for(currentNode = node; currentNode; currentNode = currentNode->next){
        if(currentNode->type == XML_ELEMENT_NODE){
            if(xmlStrcmp(currentNode->name, (const xmlChar*)"path") == 0){
                struct Path *path;
                parsePath(currentNode, &path);
                struct Geometry **geometries;
                uint16_t geometryCount = 0;
                path2Geometry(&geometries, &geometryCount, path);
                freePath(&path);

                // realloc the size for all new geometries
                uint16_t geometryOffset = svg->geometryCount;
                svg->geometryCount += geometryCount;
                svg->geometries = realloc(svg->geometries, sizeof(struct Geometry*) * svg->geometryCount);
                for(uint16_t i = 0; i < geometryCount; i++){
                    svg->geometries[i + geometryOffset] = geometries[i];
                }
            }
            // printf("%i\tnode type: element name %s\n", depth, currentNode->name);
        }

        traversElement(currentNode->children, depth + 1, svg);
    }
}
