//
//  File:       SplinesTest.cpp
//
//  Function:   Cubic spline utilities test
//
//  Copyright:  Andrew Willmott 2018
//

// Simple example that exercises the 2D api by creating a run of splines from
// a set of points, and then querying them in different ways.
//
// TODO: svg output?


#include "Splines.h"

using namespace SplineLib;

int main(int argc, const char* argv[])
{
    const Vec2f points[] =
    {
        { 1.14978, -137.651 },
        { -174.87, -78.9264 },
        { 51.2561, 90.6417 },
        { -183.824, -39.0125 },
        { 22.3669, 21.0507 },
        { -63.4809, 145.797 },
        { 16.1387, -137.128 },
        { -49.4056, 159.712 },
        { 7.183, -90.0083 },
        { -164.355, -189.983 },
        { -180.224, 42.4843 },
        { 102.902, 118.445 },
        { 8.07961, 188.87 },
        { 160.243, -126.232 },
        { -125.71, -123.763 },
        { -102.967, 18.3226 },
    };
    const int numPoints = sizeof(points) / sizeof(points[0]);
    
    cSpline2 splines[numPoints + 1];

    int numSplines = SplinesFromPoints(numPoints, points, numPoints + 1, splines);
    float sumLen = 0.0f;

    for (int i = 0; i < numSplines; i++)
    {
        float len = Length(splines[i], 0.01f);
        Bounds2f bb = ExactBounds(splines[i]);

        printf("spline %2d: length: %5.1f, bounds: [%6.1f, %6.1f], [%6.1f, %6.1f]\n", i, len, bb.min.x, bb.min.y, bb.min.x, bb.min.y);
        
        sumLen += len;
    }
    printf("\ntotal length: %g\n\n", sumLen);

    int   is = 4;
    float ts = 0.2f;
    Vec2f ps = Position (splines[is], ts);
    Vec2f vs = Velocity (splines[is], ts);
    float cs = Curvature(splines[is], ts);
    
    printf("spline %d, t = %g: position [%g, %g], velocity: [%g, %g], curvature %g \n\n", is, ts, ps.x, ps.y, vs.x, vs.y, cs);
    
    const Vec2f queryPoints[] =
    {
        { 83.7624, -122.161 },
        { 126.166, -40.4692 },
        { -36.1414, -9.16203 },
        { -61.7727, -131.043 },
    };
    
    for (Vec2f qp : queryPoints)
    {
        int index;
        float t = FindClosestPoint(qp, numSplines, splines, &index);

        Vec2f cp = Position(splines[index], t);
        
        printf("Closest point to [%6.1f, %6.1f]: index = %2d, t = %4.2f, point = [%6.1f, %6.1f]\n", qp.x, qp.y, index, t, cp.x, cp.y);
    }
    
    return 0;
}
