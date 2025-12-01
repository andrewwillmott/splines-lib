//
// SplinesTest.cpp
//
// Cubic spline utilities test
//
// Andrew Willmott
//

// Simple example that exercises the 2D/3D apis by creating a run of splines from
// a set of points, and then querying them in different ways.

#include "Splines.hpp"
#include "RotSplines.hpp"

#include <stdio.h>

using namespace SL;

void TestSplines2()
{
    printf("+ Spline2\n");

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

    Spline2 splines[numPoints + 1];

    int numSplines = SplinesFromPoints(numPoints, points, splines);
    float sumLen = 0.0f;

    printf("\ninfo:\n");
    for (int i = 0; i < numSplines; i++)
    {
        float len = Length(splines[i], 0.01f);
        Bounds2 bb = ExactBounds(splines[i]);

        printf("  %2d: length: %5.1f, bounds [%6.1f, %6.1f], [%6.1f, %6.1f]\n", i, len, bb.mMin.x, bb.mMin.y, bb.mMax.x, bb.mMax.y);

        sumLen += len;
    }
    printf("\ntotal length: %6.1f\n", sumLen);

    printf("\ninterpolation:\n");
    for (int is = 4; is < 6; is++)
        for (float ts = 0.1f; ts <= 1.0f; ts += 0.2f)
        {
            Vec2f ps = Position (splines[is], ts);
            Vec2f vs = Velocity (splines[is], ts);
            float cs = Curvature(splines[is], ts);

            printf("  i %d, t %.1f: position [%6.1f, %6.1f], velocity [%6.1f, %6.1f], curvature %6.3f\n", is, ts, ps.x, ps.y, vs.x, vs.y, cs);
        }

    const Vec2f queryPoints[] =
    {
        { 83.7624f, -122.161f },
        { 126.166f, -40.4692f },
        { -36.1414f, -9.16203f },
        { -61.7727f, -131.043f },
    };

    printf("\nclosest points:\n");
    for (Vec2f qp : queryPoints)
    {
        int index;
        float t = FindClosestPoint(qp, numSplines, splines, &index);

        Vec2f cp = Position(splines[index], t);

        printf("  to [%6.1f, %6.1f]: i %2d, t %4.2f, point [%6.1f, %6.1f]\n", qp.x, qp.y, index, t, cp.x, cp.y);
    }

    float intT[32][2];
    int   intI[32][2];

    int numIntersections = FindSplineIntersections(numSplines / 2, splines, sizeof(intT) / sizeof(intT[0]), intI, intT, 1.0f);

    printf("\nself-intersections:\n");
    for (int i = 0; i < numIntersections; i++)
        printf("  i0 %d, t0 %4.2f, i1 %d, t1 %4.2f\n", intI[i][0], intT[i][0], intI[i][1], intT[i][1]);

    printf("\n");
}

void TestSplines3()
{
    printf("+ Spline3\n");

    Vec3f points[] =
    {
        { 1.14978f , -137.651f, 11.1f },
        { -174.87f , -78.9264f, 22.2f },
        { 51.2561f , 90.6417f , 66.6f },
        { -183.824f, -39.0125f, 88.8f },
        { 22.3669f , 21.0507f , 66.6f },
        { -63.4809f, 145.797f , 44.4f },
        { 16.1387f , -137.128f, -22.2f },
        { -49.4056f, 159.712f , 11.1f },
        { 7.183f   , -90.0083f, 11.1f },
        { -164.355f, -189.983f, 33.3f },
        { -180.224f, 42.4843f , 11.1f },
        { 102.902f , 118.445f , 44.4f },
        { 8.07961f , 188.87f  , 66.6f },
        { 160.243f , -126.232f, 88.8f },
        { -125.71f , -123.763f, 55.5f },
        { -102.967f, 18.3226f , 11.1f },
    };
    const int numPoints = sizeof(points) / sizeof(points[0]);

    Spline3 splines[numPoints + 1];

    int numSplines = SplinesFromPoints(numPoints, points, splines);
    float sumLen = 0.0f;

    printf("\ninfo:\n");
    for (int i = 0; i < numSplines; i++)
    {
        float len = Length(splines[i], 0.01f);
        Bounds3 bb = ExactBounds(splines[i]);

        printf("  %2d: length %5.1f, bounds [%6.1f, %6.1f, %6.1f], [%6.1f, %6.1f, %6.1f]\n", i, len, bb.mMin.x, bb.mMin.y, bb.mMin.z, bb.mMax.x, bb.mMax.y, bb.mMax.z);

        sumLen += len;
    }
    printf("\ntotal length: %6.1f\n", sumLen);

    printf("\ninterpolation:\n");
    for (int is = 4; is < 6; is++)
        for (float ts = 0.1f; ts <= 1.0f; ts += 0.2f)
        {
            Vec3f ps = Position (splines[is], ts);
            Vec3f vs = Velocity (splines[is], ts);
            float cs = Curvature(splines[is], ts);
            float rs = Torsion  (splines[is], ts);
            printf("  i %d, t %.1f: position [%6.1f, %6.1f, %6.1f], velocity [%6.1f, %6.1f, %6.1f], curvature %6.3f, torsion %6.3f\n", is, ts, ps.x, ps.y, ps.z, vs.x, vs.y, vs.z, cs, rs);
        }

    const Vec3f queryPoints[] =
    {
        { 83.7624f , -122.161f, 22.2f },
        { 126.166f , -40.4692f, 11.1f },
        { -36.1414f, -9.16203f, 88.8f },
        { -61.7727f, -131.043f, 66.6f },
    };

    printf("\nclosest points:\n");
    for (Vec3f qp : queryPoints)
    {
        int index;
        float t = FindClosestPoint(qp, numSplines, splines, &index);

        Vec3f cp = Position(splines[index], t);

        printf("  to [%6.1f, %6.1f, %6.1f]: i %2d, t %4.2f, point [%6.1f, %6.1f, %6.1f]\n", qp.x, qp.y, qp.z, index, t, cp.x, cp.y, cp.z);
    }

    // Flatten sections so we get some self-intersections.
    points[1].z = 11.1f;
    points[2].z = 11.1f;
    points[5].z = 11.1f;
    points[6].z = 11.1f;

    SplinesFromPoints(numPoints, points, splines, 1.0f);  // tension=1 forces flat sections for intersections

    float intT[32][2];
    int   intI[32][2];

    int numIntersections = FindSplineIntersections(numSplines / 2, splines, sizeof(intT) / sizeof(intT[0]), intI, intT, 1.0f);

    printf("\nself-intersections:\n");
    for (int i = 0; i < numIntersections; i++)
        printf("  i0 %d, t0 %4.2f, i1 %d, t1 %4.2f\n", intI[i][0], intT[i][0], intI[i][1], intT[i][1]);

    printf("\n");
}

void TestRotSplines2()
{
    printf("+ RotSpline2\n");

#if 1
    const Vec2f dirs[] =
    {
        {  1,  0 },
        {  0,  1 },
        { -1,  0 },
        {  0, -1 },
        {  1,  0 },
        {  0, -1 },
    };
    const int numDirs = sizeof(dirs) / sizeof(dirs[0]);

    RotSpline2 splines[numDirs + 1];

    int numSplines = RotSplinesFromDirs(numDirs, dirs, splines);
#else
    const float angles[] =
    {
        0,
        vlf_pi * 4 + vlf_halfPi,
        vlf_pi,
        -vlf_halfPi,
        vlf_twoPi,
        vlf_twoPi - vlf_halfPi,
    };
    const int numAngles = sizeof(angles) / sizeof(angles[0]);

    RotSpline2 splines[numAngles + 1];

    int numSplines = RotSplinesFromAngles(numAngles, angles, splines);
#endif

    printf("\ninterpolation:\n");
    for (int is = 0; is < numSplines; is++)
        for (float ts = 0.1f; ts <= 1.0f; ts += 0.2f)
        {
            float w  = Rotation   (splines[is], ts);
            float wv = RotVelocity(splines[is], ts);

            Vec2f v(sinf(w), cosf(w));

            printf("  i %d, t %0.2f: direction [%6.3f, %6.3f], rot speed %.2f rad/s\n", is, ts, v.x, v.y, 2.0f * wv);
        }

    printf("\n");
}

void TestRotSplines3()
{
    printf("+ RotSpline3\n");

    const Quatf quats[] =
    {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 },
        { 0, 0, 1, 0 },
        { 0, 0, 0, 1 },
        { 0.7071, 0.7071, 0, 0 },
        { 0, 0.7071, 0.7071, 0 },
        { 0, 0, 0.7071, 0.7071 },
        { 0.7071, 0, 0, 0.7071 },
    };
    const int numQuats = sizeof(quats) / sizeof(quats[0]);

    RotSpline3 splines[numQuats + 1];

    int numSplines = RotSplinesFromQuats(numQuats, quats, splines);

    printf("\ninterpolation:\n");
    for (int is = 0; is < numSplines; is++)
        for (float ts = 0.1f; ts <= 1.0f; ts += 0.2f)
        {
            Quatf q = Rotation   (splines[is], ts);
            Vec3f w = RotVelocity(splines[is], ts);

            printf("  i %d, t %0.2f: quat [%6.3f, %6.3f, %6.3f, %6.3f], rot speed %.2f rad/s\n", is, ts, q.x, q.y, q.z, q.w, 2.0f * len(w));
        }

    printf("\n");
}

int main(int argc, const char* argv[])
{
    (void) argc; (void) argv;

    TestSplines2();
    TestSplines3();

    TestRotSplines2();
    TestRotSplines3();

    return 0;
}
