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
    // Types
    struct Vec2f { float x; float y;                    Vec2f() {}; Vec2f(float xi, float yi)                     : x(xi), y(yi)               {}; SL_VEC2F_CONVERT };
    struct Vec3f { float x; float y; float z;           Vec3f() {}; Vec3f(float xi, float yi, float zi)           : x(xi), y(yi), z(zi)        {}; SL_VEC3F_CONVERT };
    struct Vec4f { float x; float y; float z; float w;  Vec4f() {}; Vec4f(float xi, float yi, float zi, float wi) : x(xi), y(yi), z(zi), w(wi) {}; SL_VEC4F_CONVERT };
    struct Mat2f { Vec2f rows[2]; };
    struct Mat3f { Vec3f rows[3]; };
    struct Bounds2f { Vec2f min; Vec2f max; };
    struct Bounds3f { Vec3f min; Vec3f max; };

    using std::vector;


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
    cSpline2 QuadrantSpline  (const Vec2f& p, float r, int quadrant);        ///< Returns a spline representing the given quadrant (quarter circle) of radius 'r' at 'p'
    void     CircleSplines   (const Vec2f& p, float r, cSpline2 splines[4]); ///< Fills 'splines' with four splines representing a circle of radius 'r' at 'p'

    int   NumSplinesForPoints(int numPoints);     ///< Returns number of splines needed to represent the given number of points; generally n-1 except for n < 2.
    int   SplinesFromPoints  (int numPoints, const Vec2f p[], int maxSplines, cSpline2 splines[], float tension = 0.0f, size_t stride = sizeof(Vec2f));
    ///< Fills 'splines' with splines that interpolate the points in 'p', and returns the number of these splines.
    ///< 'tension' controls the interpolation -- the default value of 0 specifies Catmull-Rom splines that
    ///< guarantee tangent continuity. With +1 you get straight lines, and -1 gives more of a circular appearance.

    int   SplinesFromBezier (int numPoints, const Vec2f points[], const Vec2f hullPoints[], cSpline2 splines[], bool split = false); ///< Creates splines from the given points and Bezier hull points. If 'split' is false the splines are assumed to be continous and numPoints - 1 splinea are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.
    int   SplinesFromHermite(int numPoints, const Vec2f points[], const Vec2f tangents  [], cSpline2 splines[], bool split = false); ///< Creates splines from the given points and tangents. If 'split' is false the splines are assumed to be continous and numPoints - 1 splinea are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.

    // Queries
    Vec2f Position0(const cSpline2& spline);    ///< Starting position of spline
    Vec2f Position1(const cSpline2& spline);    ///< End position of spline

    Vec2f Velocity0(const cSpline2& spline);    ///< Starting (tangential) velocity
    Vec2f Velocity1(const cSpline2& spline);    ///< End velocity

    Vec2f Position    (const cSpline2& spline, float t); ///< Returns interpolated position
    Vec2f Velocity    (const cSpline2& spline, float t); ///< Returns interpolated velocity
    Vec2f Acceleration(const cSpline2& spline, float t); ///< Returns interpolated acceleration
    float Curvature   (const cSpline2& spline, float t); ///< Returns interpolated curvature. Curvature = 1 / r where r is the radius of the local turning circle, so 0 for flat.

    void  Frame       (const cSpline2& spline, float t, Mat2f* frame);  ///< Returns interpolated orientation matrix representing the full frame at 't'

    float LengthEstimate(const cSpline2& s, float* error);              ///< Returns estimate of length of s and optionally in 'error' the maximum error of that length.
    float Length        (const cSpline2& s, float maxError = 0.01f);    ///< Returns length of spline accurate to the given tolerance, using multiple LengthEstimate() calls.
    float Length        (const cSpline2& s, float t0, float t1, float maxError = 0.01f);    ///< Returns length of spline segment over [t0, t1].

    Bounds2f FastBounds (const cSpline2& spline);               ///< Returns fast, convervative bounds based off the convex hull of the spline.
    Bounds2f ExactBounds(const cSpline2& spline);               ///< Returns exact bounds, taking into account extrema, requires solving a quadratic

    // Subdivision
    void Split(const cSpline2& spline,          cSpline2* spline0, cSpline2* spline1);  ///< Splits 'spline' into two halves (at t = 0.5) and stores the results in 'subSplines'
    void Split(const cSpline2& spline, float t, cSpline2* spline0, cSpline2* spline1);  ///< Fills 'subSplines' with the splines corresponding to [0, t] and [t, 1]
    bool Join (const cSpline2& spline0, const cSpline2& spline1, cSpline2* spline);     ///< Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    void Split(vector<cSpline2>* splines);          ///< Subdivide each spline into two pieces
    void Split(vector<cSpline2>* splines, int n);   ///< Subdivide each spline into 'n' pieces
    void Join (vector<cSpline2>* splines);          ///< Join adjacent splines where possible

    void SubdivideForLength(vector<cSpline2>* splines, float relativeError = 0.01f);    ///< Subdivide splines to be close to linear, according to relativeError.
    void SubdivideForT     (vector<cSpline2>* splines, float error = 0.01f);            ///< Subdivide splines to be close to linear in t, i.e., arcLength

    // Nearest point
    float FindClosestPoint(const Vec2f& p, const cSpline2& spline); ///< Returns t value of the closest point on s to 'p'
    float FindClosestPoint(const Vec2f& p, int numSplines, const cSpline2 splines[], int* index); ///< Returns index of nearest spline, and 't' value of nearest point on that spline.

    struct cSubSpline2
    {
        cSpline2 mSpline;
        int      mParent;
        float    mD2;
    };

    int   FindNearbySplines(const Vec2f& p, int numSplines, const cSpline2 splines[],       vector<cSubSpline2>* nearbySplines, float* smallestFarOut = 0, int maxIter = 2);
    float FindClosestPoint (const Vec2f& p, int numSplines, const cSpline2 splines[], const vector<cSubSpline2>& nearbySplines, int* index);

    int   FindSplineIntersections(const cSpline2& spline0, const cSpline2& spline1, int maxResults, float results[][2], float tolerance = 0.1f);
    ///< Returns up to 'maxResults' intersections between the two splines, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines0, const cSpline2 splines0[], int numSplines1, const cSpline2 splines1[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    ///< Returns up to 'maxResults' intersections between the two spline lists, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines, const cSpline2 splines[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    ///< Returns up to 'maxResults' self-intersections in the given spline list, accurate to the given tolerance.

    // Linear point movement along spline set
    float AdvanceAgent(const cSpline2& spline, float t, float dl); ///< Advances 'agent' at 't' on the given spline, by dl (delta length), returning t' of the new location.

    bool  AdvanceAgent(int* index, float* t, int numSplines, const cSpline2 splines[], float dl); ///< Version of AdvanceAgent for a set of splines, assumed to be continuous. Returns false if the agent has gone off the end, in which case call ClampAgent/WrapAgent/ReverseAgent.
    void  ClampAgent  (int* index, float* t, int numSplines);    ///< Clamps agent position back to the nearest endpoint.
    void  WrapAgent   (int* index, float* t, int numSplines);    ///< Wraps the agent from the end back to the start, or vice versa.
    void  ReverseAgent(int* index, float* t);                    ///< Reverses the agent. (You must also negate the sign of dl yourself.)

    // Misc operations
    cSpline2 Reverse(const cSpline2& spline);                    ///< Reverses spline endpoints and tangents so that g(t) = f(1 - t).
    cSpline2 Offset (const cSpline2& spline, float offset);      ///< Offset spline, e.g., for stroking, +ve = to the right.

    void     Reverse(vector<cSpline2>* splines);                 ///< Reverses entire spline list
    void     Offset (vector<cSpline2>* splines, float offset);   ///< Offset splines, e.g., for stroking, +ve = to the right.



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

    int   NumSplinesForPoints(int numPoints);     ///< Returns number of splines needed to represent the given number of points; generally n-1 except for n < 2.
    int   SplinesFromPoints  (int numPoints, const Vec3f p[], int maxSplines, cSpline3 splines[], float tension = 0.0f, size_t stride = sizeof(Vec3f));
    ///< Fills 'splines' with splines that interpolate the points in 'p', and returns the number of these splines.
    ///< 'tension' controls the interpolation -- the default value of 0 specifies Catmull-Rom splines that
    ///< guarantee tangent continuity. With +1 you get straight lines, and -1 gives more of a circular appearance.

    int   SplinesFromBezier (int numPoints, const Vec3f points[], const Vec3f hullPoints[], cSpline3 splines[], bool split = false); ///< Creates splines from the given points and Bezier hull points. If 'split' is false the splines are assumed to be continous and numPoints - 1 splinea are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.
    int   SplinesFromHermite(int numPoints, const Vec3f points[], const Vec3f tangents  [], cSpline3 splines[], bool split = false); ///< Creates splines from the given points and tangents. If 'split' is false the splines are assumed to be continous and numPoints - 1 splinea are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.

    // Queries
    Vec3f Position0(const cSpline3& spline);    ///< Starting position of spline
    Vec3f Position1(const cSpline3& spline);    ///< End position of spline

    Vec3f Velocity0(const cSpline3& spline);    ///< Starting (tangential) velocity
    Vec3f Velocity1(const cSpline3& spline);    ///< End velocity

    Vec3f Position    (const cSpline3& spline, float t); ///< Returns interpolated position
    Vec3f Velocity    (const cSpline3& spline, float t); ///< Returns interpolated velocity
    Vec3f Acceleration(const cSpline3& spline, float t); ///< Returns interpolated acceleration
    float Curvature   (const cSpline3& spline, float t); ///< Returns interpolated curvature. Curvature = 1 / r where r is the radius of the local turning circle, so 0 for flat.

    void  Frame       (const cSpline3& spline, float t, Mat3f* frame); ///< Returns interpolated orientation matrix representing the full frame at 't' (z-up/y-forward/x-right)
    void  FrameSmooth (const cSpline3& spline, float t, Mat3f* frame); ///< As per Frame(), but uses the previous value of 'frame' to avoid inversions and other issues. 'frame' must be initialised on first call, to an approximation of the desired initial frame, or the identity.
    void  FrameZUp    (const cSpline3& spline, float t, Mat3f* frame); ///< Returns frame with Z-up
    void  FrameUp     (const cSpline3& spline, float t, const Vec3f& up, Mat3f* frame); ///< Returns frame with the given up vector.

    float LengthEstimate(const cSpline3& s, float* error);              ///< Returns estimate of length of s and optionally in 'error' the maximum error of that length.
    float Length        (const cSpline3& s, float maxError = 0.01f);    ///< Returns length of spline accurate to the given tolerance
    float Length        (const cSpline3& s, float t0, float t1, float maxError = 0.01f);    ///< Returns length of spline segment over [t0, t1].

    Bounds3f FastBounds (const cSpline3& spline);               ///< Returns fast, convervative bounds based off the convex hull of the spline.
    Bounds3f ExactBounds(const cSpline3& spline);               ///< Returns exact bounds, taking into account extrema, requires solving a quadratic

    // Subdivision
    void Split(const cSpline3& spline,          cSpline3* spline0, cSpline3* spline1);  ///< Splits 'spline' into two halves (at t = 0.5) and stores the results in 'subSplines'
    void Split(const cSpline3& spline, float t, cSpline3* spline0, cSpline3* spline1);  ///< Fills 'subSplines' with the splines corresponding to [0, t] and [t, 1]
    bool Join (const cSpline3& spline0, const cSpline3& spline1, cSpline3* spline);     ///< Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    void Split(vector<cSpline3>* splines);          ///< Subdivide each spline into two pieces
    void Split(vector<cSpline3>* splines, int n);   ///< Subdivide each spline into 'n' pieces
    void Join (vector<cSpline3>* splines);          ///< Join adjacent splines where possible

    void SubdivideForLength(vector<cSpline3>* splines, float relativeError = 0.01f);    ///< Subdivide splines to be close to linear, according to relativeError.
    void SubdivideForT     (vector<cSpline3>* splines, float error = 0.01f);            ///< Subdivide splines to be close to linear in t, i.e., arcLength

    // Nearest point
    float FindClosestPoint(const Vec3f& p, const cSpline3& spline); ///< Returns t value of the closest point on s to 'p'
    float FindClosestPoint(const Vec3f& p, int numSplines, const cSpline3 splines[], int* index); ///< Returns index of nearest spline, and 't' value of nearest point on that spline.

    struct cSubSpline3
    {
        cSpline3 mSpline;
        int      mParent;
        float    mD2;
    };

    int   FindNearbySplines(const Vec3f& p, int numSplines, const cSpline3 splines[],       vector<cSubSpline3>* nearbySplines, float* smallestFarOut = 0, int maxIter = 2);
    float FindClosestPoint (const Vec3f& p, int numSplines, const cSpline3 splines[], const vector<cSubSpline3>& nearbySplines, int* index);

    int   FindSplineIntersections(const cSpline3& spline0, const cSpline3& spline1, int maxResults, float results[][2], float tolerance = 0.1f);
    ///< Returns up to 'maxResults' intersections between the two splines, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines0, const cSpline3 splines0[], int numSplines1, const cSpline3 splines1[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    ///< Returns up to 'maxResults' intersections between the two spline lists, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines, const cSpline3 splines[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    ///< Returns up to 'maxResults' self-intersections in the given spline list, accurate to the given tolerance.

    // Linear point movement along spline set
    float AdvanceAgent(const cSpline3& spline, float t, float dl); ///< Advances 'agent' at 't' on the given spline, by dl (delta length), returning t' of the new location.

    bool  AdvanceAgent(int* index, float* t, int numSplines, const cSpline3 splines[], float dl); ///< Version of AdvanceAgent for a set of splines, assumed to be continuous. Returns false if the agent has gone off the end, in which case call ClampAgent/WrapAgent/ReverseAgent.
//  void  ClampAgent  (int* index, float* t, int numSplines);    ///< Clamps agent position back to the nearest endpoint.
//  void  WrapAgent   (int* index, float* t, int numSplines);    ///< Wraps the agent from the end back to the start, or vice versa.
//  void  ReverseAgent(int* index, float* t);                    ///< Reverses the agent. (You must also negate the sign of dl yourself.)

    // Misc operations
    cSpline3 Reverse(const cSpline3& spline);                    ///< Reverses spline endpoints and tangents so that g(t) = f(1 - t).
    cSpline3 Offset (const cSpline3& spline, float offset);      ///< Offset spline, e.g., for stroking, +ve = to the right.

    void     Reverse(vector<cSpline3>* splines);                 ///< Reverses entire spline list
    void     Offset (vector<cSpline3>* splines, float offset);   ///< Offset splines, e.g., for stroking, +ve = to the right.
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
