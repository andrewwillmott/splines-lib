//
//  File:       SplinesTest.cpp
//
//  Function:   Cubic spline utilities test
//
//  Copyright:  Andrew Willmott 2018
//

// Simple example that exercises the 2D api by creating a run of splines from
// a set of points, and then querying them in different ways.


#include "Splines.h"

#include <stdio.h>

using namespace SplineLib;

int main(int argc, const char* argv[])
{
    const Vec2f points[] =
    {
        { 1.14978f, -137.651f },
        { -174.87f, -78.9264f },
        { 51.2561f, 90.6417f },
        { -183.824f, -39.0125f },
        { 22.3669f, 21.0507f },
        { -63.4809f, 145.797f },
        { 16.1387f, -137.128f },
        { -49.4056f, 159.712f },
        { 7.183f, -90.0083f },
        { -164.355f, -189.983f },
        { -180.224f, 42.4843f },
        { 102.902f, 118.445f },
        { 8.07961f, 188.87f },
        { 160.243f, -126.232f },
        { -125.71f, -123.763f },
        { -102.967f, 18.3226f },
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
        { 83.7624f, -122.161f },
        { 126.166f, -40.4692f },
        { -36.1414f, -9.16203f },
        { -61.7727f, -131.043f },
    };

    for (Vec2f qp : queryPoints)
    {
        int index;
        float t = FindClosestPoint(qp, numSplines, splines, &index);

        Vec2f cp = Position(splines[index], t);
        
        printf("Closest point to [%6.1f, %6.1f]: index = %2d, t = %4.2f, point = [%6.1f, %6.1f]\n", qp.x, qp.y, index, t, cp.x, cp.y);
    }

    float intT[32][2];
    int   intI[32][2];

    int numIntersections = FindSplineIntersections(numSplines / 2, splines, sizeof(intT) / sizeof(intT[0]), intI, intT, 1.0f);

    printf("\nSelf-intersections:\n");
    for (int i = 0; i < numIntersections; i++)
        printf("    i0 = %2d, t0 = %4.2f, i1 = %2d, t1 = %4.2f\n", intI[i][0], intT[i][0], intI[i][1], intT[i][1]);

    return 0;
}
