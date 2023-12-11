#include "mathFunc.h"

#include "math.h"

double rad(double r){
    return (r * 180.0) / ((double)M_PI);
}

double angleCalc(double x1, double y1, double x2, double y2){
    double innerProduct = x1 * x2 + y1 * y2;
    double lengthVectro1 = hypot(x1, y1);
    double lengthVectro2 = hypot(x2, y2);
    double r = acos(innerProduct / (lengthVectro1 * lengthVectro2));

    if(x1 * y2 - y1 * x2 < 0.0){
        r *= -1.0;
    }

    return r;
}

void normalizeAngle360(double *angle){
    while(*angle < 0.0){
        *angle += 360.0;
    }
    while(360.0 <= *angle){
        *angle -= 360.0;
    }
}

void rotate(double *x, double *y, double theta){
    double xPrime = *x * cos(theta) - *y * sin(theta);
    double yPrime = *x * sin(theta) + *y * cos(theta);
    *x = xPrime;
    *y = yPrime;
}
