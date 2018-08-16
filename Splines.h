//
//  File:       Splines.h
//
//  Function:   Various cubic spline utilities
//
//  Copyright:  Andrew Willmott 2018
//

#ifndef SPLINES_H
#define SPLINES_H

#include <math.h>
#include <vector>

#ifndef SL_ASSERT
    #define SL_ASSERT(X)
#endif

#ifndef SL_VEC2F_CONVERT
    #define SL_VEC2F_CONVERT    // E.g, #define SL_VEC2F_CONVERT Vec2f(MyV2 v2) : x(v2.x), y(v2.y) {}; operator MyV2() const { return { x, y }; } 
    #define SL_VEC3F_CONVERT
    #define SL_VEC4F_CONVERT
#endif

namespace SplineLib
{
    struct Vec2f { float x; float y;                    Vec2f() {}; Vec2f(float xi, float yi)                     : x(xi), y(yi)               {}; SL_VEC2F_CONVERT };
    struct Vec3f { float x; float y; float z;           Vec3f() {}; Vec3f(float xi, float yi, float zi)           : x(xi), y(yi), z(zi)        {}; SL_VEC3F_CONVERT };
    struct Vec4f { float x; float y; float z; float w;  Vec4f() {}; Vec4f(float xi, float yi, float zi, float wi) : x(xi), y(yi), z(zi), w(wi) {}; SL_VEC4F_CONVERT };

    struct Bounds2f { Vec2f min; Vec2f max; };
    struct Bounds3f { Vec3f min; Vec3f max; };


    // 2D
    struct cSpline2
    {
        Vec4f xb;   // x cubic bezier coefficients
        Vec4f yb;   // y cubic bezier coefficients
    };

    // Spline creation
    cSpline2 BezierSpline    (const Vec2f& p0, const Vec2f& p1, const Vec2f& p2, const Vec2f& p3);  ///< Return Bezier spline from p0 to p3 with guide points p1, p2
    cSpline2 HermiteSpline   (const Vec2f& p0, const Vec2f& p1, const Vec2f& v0, const Vec2f& v1);  ///< Return Hermite spline from p0 to p1 with corresponding tangents v0, v1.
    cSpline2 CatmullRomSpline(const Vec2f& p0, const Vec2f& p1, const Vec2f& p2, const Vec2f& p3);  ///< Return Catmull-Rom spline passing through p1 and p2, with tangents affected by p0 and p3.
    
    int NumSplinesForPoints(int numPoints);     ///< Returns number of splines needed to represent the given number of points; generally n-1 except for n < 2.
    int SplinesFromPoints  (int numPoints, const Vec2f p[], int maxSplines, cSpline2 splines[], float tension = 0.0f, size_t stride = sizeof(Vec2f));
    ///< Fills 'splines' with splines that interpolate the points in 'p', and returns the number of these splines.
    ///< 'tension' controls the interpolation -- the default value of 0 specifies Catmull-Rom splines that
    ///< guarantee tangent continuity. With +1 you get straight lines, and -1 gives more of a circular appearance.

    // Queries
    Vec2f Position0(const cSpline2& spline);    ///< Starting position of spline
    Vec2f Position1(const cSpline2& spline);    ///< End position of spline
    
    Vec2f Velocity0(const cSpline2& spline);    ///< Starting (tangential) velocity
    Vec2f Velocity1(const cSpline2& spline);    ///< End velocity
    
    Vec2f Position    (const cSpline2& spline, float t); ///< Returns interpolated position
    Vec2f Velocity    (const cSpline2& spline, float t); ///< Returns interpolated velocity
    Vec2f Acceleration(const cSpline2& spline, float t); ///< Returns interpolated acceleration
    float Curvature   (const cSpline2& spline, float t); ///< Returns interpolated curvature. Curvature = 1 / r where r is the radius of the local turning circle, so 0 for flat.
    
    float LengthEstimate(const cSpline2& s, float* error);      ///< Returns estimate of length of s and optionally in 'error' the maximum error of that length.
    float Length        (const cSpline2& s, float maxError);    ///< Returns length of spline accurate to the given tolerance, using multiple LengthEstimate() calls.

    Bounds2f FastBounds (const cSpline2& spline);               ///< Returns fast, convervative bounds based off the convex hull of the spline.
    Bounds2f ExactBounds(const cSpline2& spline);               ///< Returns exact bounds, taking into account extrema, requires solving a quadratic

    // Subdivision
    void Split(const cSpline2& spline,          cSpline2* spline0, cSpline2* spline1);  ///< Splits 'spline' into two halves (at t = 0.5) and stores the results in 'subSplines'
    void Split(const cSpline2& spline, float t, cSpline2* spline0, cSpline2* spline1);  ///< Fills 'subSplines' with the splines corresponding to [0, t] and [t, 1]
    bool Join (const cSpline2& spline0, const cSpline2& spline1, cSpline2* spline);     ///< Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    // Nearest point
    float FindClosestPoint(const Vec2f& p, const cSpline2& spline); ///< Returns t value of the closest point on s to 'p'
    float FindClosestPoint(const Vec2f& p, int numSplines, const cSpline2 splines[], int* index); ///< Returns index of nearest spline, and 't' value of nearest point on that spline.

    struct cSubSpline2
    {
        cSpline2 mSpline;
        float    mD2;
        int      mParent;
    };

    int   FindNearbySplines(const Vec2f& p, int numSplines, const cSpline2 splines[],       std::vector<cSubSpline2>* nearbySplines, float* smallestFarOut = 0, int maxIter = 2);
    float FindClosestPoint (const Vec2f& p, int numSplines, const cSpline2 splines[], const std::vector<cSubSpline2>& nearbySplines, int* index);

    // Linear point movement along spline set
    float AdvanceAgent(const cSpline2& spline, float t, float dl);  ///< Advances 'agent' at 't' on the given spline, by dl (delta length), returning t' of the new location.
    
    bool AdvanceAgent(int* index, float* t, int numSplines, const cSpline2 splines[], float dl); ///< Version of AdvanceAgent for a set of splines, assumed to be continuous. Returns false if the agent has gone off the end, in which case call ClampAgent/WrapAgent/ReverseAgent.
    void ClampAgent  (int* index, float* t, int numSplines);    ///< Clamps agent position back to the nearest endpoint.
    void WrapAgent   (int* index, float* t, int numSplines);    ///< Wraps the agent from the end back to the start, or vice versa.
    void ReverseAgent(int* index, float* t);                    ///< Reverses the agent. (You must also negate the sign of dl yourself.)

    
    // 3D
    struct cSpline3
    {
        Vec4f xb;   // x cubic bezier coefficients
        Vec4f yb;   // y cubic bezier coefficients
        Vec4f zb;   // z cubic bezier coefficients
    };

    // Spline creation
    cSpline3 BezierSpline    (const Vec3f& p0, const Vec3f& p1, const Vec3f& p2, const Vec3f& p3);  ///< Return Bezier spline from p0 to p3 with guide points p1, p2
    cSpline3 HermiteSpline   (const Vec3f& p0, const Vec3f& p1, const Vec3f& v0, const Vec3f& v1);  ///< Return Hermite spline from p0 to p1 with corresponding tangents v0, v1.
    cSpline3 CatmullRomSpline(const Vec3f& p0, const Vec3f& p1, const Vec3f& p2, const Vec3f& p3);  ///< Return Catmull-Rom spline passing through p1 and p2, with tangents affected by p0 and p3.

    int NumSplinesForPoints(int numPoints);     ///< Returns number of splines needed to represent the given number of points; generally n-1 except for n < 2.
    int SplinesFromPoints  (int numPoints, const Vec3f p[], int maxSplines, cSpline3 splines[], float tension = 0.0f, size_t stride = sizeof(Vec3f));
    ///< Fills 'splines' with splines that interpolate the points in 'p', and returns the number of these splines.
    ///< 'tension' controls the interpolation -- the default value of 0 specifies Catmull-Rom splines that
    ///< guarantee tangent continuity. With +1 you get straight lines, and -1 gives more of a circular appearance.

    // Queries
    Vec3f Position0(const cSpline3& spline);    ///< Starting position of spline
    Vec3f Position1(const cSpline3& spline);    ///< End position of spline

    Vec3f Velocity0(const cSpline3& spline);    ///< Starting (tangential) velocity
    Vec3f Velocity1(const cSpline3& spline);    ///< End velocity

    Vec3f Position    (const cSpline3& spline, float t); ///< Returns interpolated position
    Vec3f Velocity    (const cSpline3& spline, float t); ///< Returns interpolated velocity
    Vec3f Acceleration(const cSpline3& spline, float t); ///< Returns interpolated acceleration
    float Curvature   (const cSpline3& spline, float t); ///< Returns interpolated curvature. Curvature = 1 / r where r is the radius of the local turning circle, so 0 for flat.

    float LengthEstimate(const cSpline3& s, float* error);      ///< Returns estimate of length of s and optionally in 'error' the maximum error of that length.
    float Length        (const cSpline3& s, float maxError);    ///< Returns length of spline accurate to the given tolerance

    Bounds3f FastBounds (const cSpline3& spline);               ///< Returns fast, convervative bounds based off the convex hull of the spline.
    Bounds3f ExactBounds(const cSpline3& spline);               ///< Returns exact bounds, taking into account extrema, requires solving a quadratic

    // Subdivision
    void Split(const cSpline3& spline,          cSpline3* spline0, cSpline3* spline1);  ///< Splits 'spline' into two halves (at t = 0.5) and stores the results in 'subSplines'
    void Split(const cSpline3& spline, float t, cSpline3* spline0, cSpline3* spline1);  ///< Fills 'subSplines' with the splines corresponding to [0, t] and [t, 1]
    bool Join (const cSpline3& spline0, const cSpline3& spline1, cSpline3* spline);     ///< Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    // Nearest point
    float FindClosestPoint(const cSpline3& spline, const Vec3f& p); ///< Returns t value of the closest point on s to 'p'
    float FindClosestPoint(const Vec3f& p, int numSplines, const cSpline3 splines[], int* index); ///< Returns index of nearest spline, and 't' value of nearest point on that spline.

    struct cSubSpline3
    {
        cSpline3 mSpline;
        float    mD2;
        int      mParent;
    };

    int   FindNearbySplines(const Vec3f& p, int numSplines, const cSpline3 splines[],       std::vector<cSubSpline3>* nearbySplines, float* smallestFarOut = 0, int maxIter = 2);
    float FindClosestPoint (const Vec3f& p, int numSplines, const cSpline3 splines[], const std::vector<cSubSpline3>& nearbySplines, int* index);

    // Linear point movement along spline set
    float AdvanceAgent(const cSpline3& spline, float t, float dl);  ///< Advances 'agent' at 't' on the given spline, by dl (delta length), returning t' of the new location.
    
    bool AdvanceAgent(int* index, float* t, int numSplines, const cSpline3 splines[], float dl); ///< Version of AdvanceAgent for a set of splines, assumed to be continuous. Returns false if the agent has gone off the end, in which case call ClampAgent/WrapAgent/ReverseAgent.
//  void ClampAgent  (int* index, float* t, int numSplines);    ///< Clamps agent position back to the nearest endpoint.
//  void WrapAgent   (int* index, float* t, int numSplines);    ///< Wraps the agent from the end back to the start, or vice versa.
//  void ReverseAgent(int* index, float* t);                    ///< Reverses the agent. (You must also negate the sign of dl yourself.)
}


// --- Inlines -----------------------------------------------------------------

inline int SplineLib::NumSplinesForPoints(int numPoints)
{
    if (numPoints < 2)
        return numPoints;

    return numPoints - 1;
}

// 2D

inline SplineLib::Vec2f SplineLib::Position0(const cSpline2& spline)
{
    return Vec2f(spline.xb.x, spline.yb.x);
}
inline SplineLib::Vec2f SplineLib::Position1(const cSpline2& spline)
{
    return Vec2f(spline.xb.w, spline.yb.w);
}

inline SplineLib::Vec2f SplineLib::Velocity0(const cSpline2& spline)
{
    return Vec2f
    (
        3.0f * (spline.xb.y - spline.xb.x),
        3.0f * (spline.yb.y - spline.yb.x)
    );
}
inline SplineLib::Vec2f SplineLib::Velocity1(const cSpline2& spline)
{
    return Vec2f
    (
        3.0f * (spline.xb.w - spline.xb.z),
        3.0f * (spline.yb.w - spline.yb.z)
    );
}

inline void SplineLib::ClampAgent(int* index, float* t, int numSplines)
{
    if (*index < 0)
    {
        *index = 0;
        *t = 0.0f;
    }
    else if (*index >= numSplines)
    {
        *index = numSplines - 1;
        *t = 1.0f;
    } 
    else if (*t < 0.0f)
        *t = 0.0f;
    else if (*t > 1.0f)
        *t = 1.0f;
}

inline void SplineLib::WrapAgent(int* indexInOut, float* tInOut, int numSplines)
{
    int& index = *indexInOut;
    float& t = *tInOut;

    SL_ASSERT(!IsNAN(t));
    SL_ASSERT(index == 0 || index == numSplines - 1);

    t -= floorf(t);
    index ^= numSplines - 1;
}

inline void SplineLib::ReverseAgent(int* , float* t)
{
    *t = ceilf(*t) - *t;
}

// 3D

inline SplineLib::Vec3f SplineLib::Position0(const cSpline3& spline)
{
    return Vec3f(spline.xb.x, spline.yb.x, spline.zb.x);
}
inline SplineLib::Vec3f SplineLib::Position1(const cSpline3& spline)
{
    return Vec3f(spline.xb.w, spline.yb.w, spline.zb.w);
}

inline SplineLib::Vec3f SplineLib::Velocity0(const cSpline3& spline)
{
    return Vec3f
    (
        3.0f * (spline.xb.y - spline.xb.x),
        3.0f * (spline.yb.y - spline.yb.x),
        3.0f * (spline.zb.y - spline.zb.x)
    );
}
inline SplineLib::Vec3f SplineLib::Velocity1(const cSpline3& spline)
{
    return Vec3f
    (
        3.0f * (spline.xb.w - spline.xb.z),
        3.0f * (spline.yb.w - spline.yb.z),
        3.0f * (spline.zb.w - spline.zb.z)
    );
}

#endif
