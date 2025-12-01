//
// Splines.hpp
//
// Cubic spline utilities for 1D/2D/3D
//
// Andrew Willmott
//

#ifndef SL_SPLINES_H
#define SL_SPLINES_H

#include "VLMini.hpp"
#include <vector>

#ifndef SL_ASSERT
    #define SL_ASSERT(X)
#endif

namespace SL
{
    typedef Vec2f Bounds1;
    struct Bounds2 { Vec2f mMin; Vec2f mMax; };
    struct Bounds3 { Vec3f mMin; Vec3f mMax; };

    // 1D
    typedef Vec4f Spline1;

    // Spline creation
    Spline1 BezierSpline       (float p0, float p1, float p2, float p3); // Returns Bezier spline from p0 to p3 with guide points p1, p2
    Spline1 HermiteSpline      (float p0, float p1, float v0, float v1); // Returns Hermite spline from p0 to p1 with corresponding tangents (slopes) v0, v1.
    Spline1 CatmullRomSpline   (float p0, float p1, float p2, float p3); // Returns Catmull-Rom spline passing through p1 and p2, with tangents affected by p0 and p3.
    Spline1 InterpolatingSpline(float p0, float p1, float p2, float p3); // Returns spline that interpolates p0, p1, p2, p3 at t=0,1/3,2/3,1
    Spline1 LineSpline         (float p0, float p1);                     // Returns a spline representing the line p0_p1
    Spline1 CubicSpline        (const Vec4f& abcd);                      // Returns a spline representing a + bx + cx^2 + dx^3

    int   NumSplinesForPoints(int numPoints);     // Returns number of splines needed to represent the given number of points; generally n-1 except for n < 2.
    int   SplinesFromPoints  (int numPoints, const float p[], Spline1 splines[], float tension = 0.0f, size_t stride = sizeof(float));
    // Fills 'splines' with splines that interpolate the values in 'p', and returns the number of these splines.
    // 'tension' controls the interpolation -- the default value of 0 specifies Catmull-Rom splines that guarantee tangent continuity. With +1 you get linear interpolation, and -1 gives exaggerated smoothing.

    int   SplinesFromPointsDynamic(int numPoints, const float p[], Spline1 splines[], float tension = 0.0f, float overshootRatio = 0.5f, size_t stride = sizeof(Vec2f));
    // A version of SplinesFromPoints that adjusts tension per point to avoid overshoots/loops
    // caused by disparate segment sizes.

    int   SplinesFromBezier (int numPoints, const float points[], const float hullPoints[], Spline1 splines[], bool split = false); // Creates splines from the given points and Bezier hull points. If 'split' is false the splines are assumed to be continuous and numPoints - 1 splines are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.
    int   SplinesFromHermite(int numPoints, const float points[], const float tangents  [], Spline1 splines[], bool split = false); // Creates splines from the given points and tangents. If 'split' is false the splines are assumed to be continuous and numPoints - 1 splines are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.

    void  MakeMonotonic(int numSplines, Spline1 splines[], bool closed = false);
    // Makes splines monotonic (internal spline values always within the end points) by adjusting tangents. Useful for situations where overshoots lead to invalid values.

    // Queries
    float Position0(const Spline1& spline); // Start of spline
    float Position1(const Spline1& spline); // End of spline
    float Velocity0(const Spline1& spline); // Start velocity (dx/dt)
    float Velocity1(const Spline1& spline); // End velocity

    float Position    (const Spline1& spline, float t); // Returns position at 't' (x)
    float Velocity    (const Spline1& spline, float t); // Returns velocity at 't' (x')
    float Acceleration(const Spline1& spline, float t); // Returns acceleration at 't' (x'')
    float Jerk        (const Spline1& spline, float t); // Returns jerk at 't' (x''') -- for cubic splines this is constant

    Bounds1 FastBounds (const Spline1& spline); // Returns fast, convervative bounds based off the convex hull of the spline.
    Bounds1 ExactBounds(const Spline1& spline); // Returns exact bounds, taking into account extrema, requires solving a quadratic

    // Conversion
    Vec4f CubicCoeffs(const Spline1& spline);   // Returns the spline as cubic coefficients, [a, b, c, d]

    // Subdivision
    void Split(const Spline1& spline,          Spline1* spline0, Spline1* spline1); // Splits 'spline' into two halves (at t = 0.5) and stores the results in 'subSplines'
    void Split(const Spline1& spline, float t, Spline1* spline0, Spline1* spline1); // Fills 'subSplines' with the splines corresponding to [0, t] and [t, 1]
    bool Join (const Spline1& spline0, const Spline1& spline1, Spline1* spline);    // Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    Spline1 Trim(const Spline1& spline, float t0, float t1); // Returns subsection of 'spline' between t0 and t1


    // 2D
    struct Spline2
    {
        Spline1 xb;   // x cubic bezier coefficients
        Spline1 yb;   // y cubic bezier coefficients
    };

    // Spline creation
    Spline2 BezierSpline       (Vec2f p0, Vec2f p1, Vec2f p2, Vec2f p3); // Returns Bezier spline from p0 to p3 with guide points p1, p2
    Spline2 HermiteSpline      (Vec2f p0, Vec2f p1, Vec2f v0, Vec2f v1); // Returns Hermite spline from p0 to p1 with corresponding tangents v0, v1.
    Spline2 CatmullRomSpline   (Vec2f p0, Vec2f p1, Vec2f p2, Vec2f p3); // Returns Catmull-Rom spline passing through p1 and p2, with tangents affected by p0 and p3.
    Spline2 InterpolatingSpline(Vec2f p0, Vec2f p1, Vec2f p2, Vec2f p3); // Returns spline that interpolates p0, p1, p2, p3

    Spline2 LineSpline      (Vec2f p0, Vec2f p1);                   // Returns a spline representing the line p0_p1
    Spline2 QuadrantSpline  (Vec2f p, float r, int quadrant);       // Returns a spline representing the given quadrant (quarter circle) of radius 'r' at 'p'
    void    CircleSplines   (Vec2f p, float r, Spline2 splines[4]); // Fills 'splines' with four splines representing a circle of radius 'r' at 'p'

    int   NumSplinesForPoints(int numPoints);     // Returns number of splines needed to represent the given number of points; generally n-1 except for n < 2.
    int   SplinesFromPoints  (int numPoints, const Vec2f p[], Spline2 splines[], float tension = 0.0f, size_t stride = sizeof(Vec2f));
    // Fills 'splines' with splines that interpolate the points in 'p', and returns the number of these splines.
    // 'tension' controls the interpolation -- the default value of 0 specifies Catmull-Rom splines that
    // guarantee tangent continuity. With +1 you get straight lines, and -1 gives more of a circular appearance.

    int   SplinesFromPointsDynamic(int numPoints, const Vec2f p[], Spline2 splines[], float tension = 0.0f, float overshootRatio = 0.5f, size_t stride = sizeof(Vec2f));
    // A version of SplinesFromPoints that adjusts tension per point to avoid overshoots/loops
    // caused by disparate segment sizes.

    int   SplinesFromBezier (int numPoints, const Vec2f points[], const Vec2f hullPoints[], Spline2 splines[], bool split = false); // Creates splines from the given points and Bezier hull points. If 'split' is false the splines are assumed to be continuous and numPoints - 1 splines are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.
    int   SplinesFromHermite(int numPoints, const Vec2f points[], const Vec2f tangents  [], Spline2 splines[], bool split = false); // Creates splines from the given points and tangents. If 'split' is false the splines are assumed to be continuous and numPoints - 1 splines are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.

    void  MakeMonotonic(int numSplines, Spline2 splines[], bool closed = false);
    // Makes splines monotonic (internal spline values always within the end points) by adjusting tangents. Useful for situations where overshoots lead to invalid values.

    // Queries
    Vec2f Position0(const Spline2& spline); // Start of spline
    Vec2f Position1(const Spline2& spline); // End of spline
    Vec2f Velocity0(const Spline2& spline); // Start velocity (tangent x dp/dt)
    Vec2f Velocity1(const Spline2& spline); // End velocity

    Vec2f Position    (const Spline2& spline, float t); // Returns position at 't'
    Vec2f Velocity    (const Spline2& spline, float t); // Returns velocity at 't'
    Vec2f Acceleration(const Spline2& spline, float t); // Returns acceleration at 't'
    Vec2f Jerk        (const Spline2& spline, float t); // Returns jerk at 't' -- for cubic splines this is constant

    float Curvature   (const Spline2& spline, float t); // Returns curvature at 't'. Curvature = 1 / r, where r is the radius of the osculating circle, so 0 on a straight path.

    Mat2f Frame       (const Spline2& spline, float t);               // Returns the Frenet frame at 't': x=forward, y=left/right, depending on which way the path is curving.
    void  Frame       (const Spline2& spline, float t, Mat2f* frame); // Updates 'frame' using the Frenet frame at 't', maintaining consistency by avoiding any axis flips.
    Mat2f FrameX      (const Spline2& spline, float t);               // Returns frame with x forwards and y=left

    float LengthEstimate(const Spline2& s, float* error);             // Returns estimate of length of s and optionally in 'error' the maximum error of that length.
    float Length        (const Spline2& s, float maxError = 0.01f);   // Returns length of spline accurate to the given tolerance, using multiple LengthEstimate() calls.
    float Length        (const Spline2& s, float t0, float t1, float maxError = 0.01f); // Returns length of spline segment over [t0, t1].

    Bounds2 FastBounds (const Spline2& spline); // Returns fast, convervative bounds based off the convex hull of the spline.
    Bounds2 ExactBounds(const Spline2& spline); // Returns exact bounds, taking into account extrema, requires solving a quadratic

    // Line conversion
    constexpr int kMaxEmitLines = 128;
    typedef void tEmitLinesFunc2(int numLines, Vec2f lines[][2], void* context); // Called with up to kMaxEmitLines lines at a time
    int SplinesToLinesLinear  (int numSplines, const Spline2 splines[], tEmitLinesFunc2 emitLines, void* context, float step = 0.05f); // Convert given splines into lines using even steps in 't'
    int SplinesToLinesAdaptive(int numSplines, const Spline2 splines[], tEmitLinesFunc2 emitLines, void* context, float tolerance = 0.05f); // Convert given splines into lines adaptively, with given tolerance in world units.

    typedef void tEmitLinesExFunc2(int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* context); // Called with up to kMaxEmitLines lines at a time, plus t values and owning spline.
    int SplinesToLinesLinear  (int numSplines, const Spline2 splines[], tEmitLinesExFunc2 emitLines, void* context, float step = 0.05f);
    int SplinesToLinesAdaptive(int numSplines, const Spline2 splines[], tEmitLinesExFunc2 emitLines, void* context, float tolerance = 0.05f);

    // Subdivision
    void Split(const Spline2& spline,          Spline2* spline0, Spline2* spline1); // Splits 'spline' into two halves (at t = 0.5) and stores the results in 'subSplines'
    void Split(const Spline2& spline, float t, Spline2* spline0, Spline2* spline1); // Fills 'subSplines' with the splines corresponding to [0, t] and [t, 1]
    bool Join (const Spline2& spline0, const Spline2& spline1, Spline2* spline);    // Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    void Split(std::vector<Spline2>* splines);          // Subdivide each spline into two pieces
    void Split(std::vector<Spline2>* splines, int n);   // Subdivide each spline into 'n' pieces
    void Join (std::vector<Spline2>* splines);          // Join adjacent splines where possible

    Spline2 Trim(const Spline2& spline, float t0, float t1); // Returns subsection of 'spline' between t0 and t1

    void SubdivideForLength(std::vector<Spline2>* splines, float relativeError = 0.01f); // Subdivide splines to be close to linear, according to relativeError.
    void SubdivideForT     (std::vector<Spline2>* splines, float error = 0.01f);         // Subdivide splines to be close to linear in t, i.e., arcLength

    void SubdivideForLengthRatio(std::vector<Vec2f>& positions, float maxRatio = 8.0f);  // Insert additional positions as necessary to ensure no segment is more than 'maxRatio' times its neighbours. Useful as preconditioner for SplinesFromPoints.

    // Nearest point
    float FindClosestPoint(Vec2f p, const Spline2& spline); // Returns t value of the closest point on s to 'p'
    float FindClosestPoint(Vec2f p, int numSplines, const Spline2 splines[], int* index); // Returns index of nearest spline, and 't' value of nearest point on that spline.

    struct SubSpline2
    {
        Spline2 mSpline;
        int     mParent;
        float   mD2;
    };

    int   FindNearbySplines(Vec2f p, int numSplines, const Spline2 splines[],       std::vector<SubSpline2>* nearbySplines, float* smallestFarOut = 0, int maxIter = 2);
    float FindClosestPoint (Vec2f p, int numSplines, const Spline2 splines[], const std::vector<SubSpline2>& nearbySplines, int* index);

    int   FindSplineIntersections(const Spline2& spline0, const Spline2& spline1, int maxResults, float results[][2], float tolerance = 0.1f);
    // Returns up to 'maxResults' intersections between the two splines, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines0, const Spline2 splines0[], int numSplines1, const Spline2 splines1[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    // Returns up to 'maxResults' intersections between the two spline lists, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines, const Spline2 splines[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    // Returns up to 'maxResults' self-intersections in the given spline list, accurate to the given tolerance.

    // Linear point movement along spline set
    float AdvanceAgent(const Spline2& spline, float t, float dl); // Advances 'agent' at 't' on the given spline, by dl (delta length), returning t' of the new location.

    bool  AdvanceAgent(int* index, float* t, int numSplines, const Spline2 splines[], float dl); // Version of AdvanceAgent for a set of splines, assumed to be continuous. Returns false if the agent has gone off the end, in which case call ClampAgent/WrapAgent/ReverseAgent.
    void  ClampAgent  (int* index, float* t, int numSplines);    // Clamps agent position back to the nearest endpoint.
    void  WrapAgent   (int* index, float* t, int numSplines);    // Wraps the agent from the end back to the start, or vice versa.
    void  ReverseAgent(int* index, float* t);                    // Reverses the agent. (You must also negate the sign of dl yourself.)

    // Misc operations
    Spline2 Reverse(const Spline2& spline);                      // Reverses spline endpoints and tangents so that g(t) = f(1 - t).
    Spline2 Offset (const Spline2& spline, float offset);        // Offset spline, e.g., for stroking, +ve = to the right.

    void  Reverse(std::vector<Spline2>* splines);               // Reverses entire spline list
    void  Offset (std::vector<Spline2>* splines, float offset); // Offset splines, e.g., for stroking, +ve = to the right.


    // 3D
    struct Spline3
    {
        Spline1 xb;   // x cubic bezier coefficients
        Spline1 yb;   // y cubic bezier coefficients
        Spline1 zb;   // z cubic bezier coefficients
    };

    // Spline creation
    Spline3 BezierSpline       (Vec3f p0, Vec3f p1, Vec3f p2, Vec3f p3); // Returns Bezier spline from p0 to p3 with guide points p1, p2
    Spline3 HermiteSpline      (Vec3f p0, Vec3f p1, Vec3f v0, Vec3f v1); // Returns Hermite spline from p0 to p1 with corresponding tangents v0, v1.
    Spline3 CatmullRomSpline   (Vec3f p0, Vec3f p1, Vec3f p2, Vec3f p3); // Returns Catmull-Rom spline passing through p1 and p2, with tangents affected by p0 and p3.
    Spline3 InterpolatingSpline(Vec3f p0, Vec3f p1, Vec3f p2, Vec3f p3); // Returns spline that interpolates p0, p1, p2, p3.

    Spline3 LineSpline      (Vec3f p0, Vec3f p1);                   // Returns a spline representing the line p0_p1
    Spline3 QuadrantSpline  (Vec3f p, float r, int quadrant);       // Returns a spline representing the given quadrant (quarter circle) of radius 'r' at 'p'
    void    CircleSplines   (Vec3f p, float r, Spline3 splines[4]); // Fills 'splines' with four splines representing a circle of radius 'r' at 'p'

    int   NumSplinesForPoints(int numPoints);     // Returns number of splines needed to represent the given number of points; generally n-1 except for n < 2.
    int   SplinesFromPoints  (int numPoints, const Vec3f p[], Spline3 splines[], float tension = 0.0f, size_t stride = sizeof(Vec3f));
    // Fills 'splines' with splines that interpolate the points in 'p', and returns the number of these splines.
    // 'tension' controls the interpolation -- the default value of 0 specifies Catmull-Rom splines that
    // guarantee tangent continuity. With +1 you get straight lines, and -1 gives more of a circular appearance.

    int   SplinesFromPointsDynamic(int numPoints, const Vec3f p[], Spline3 splines[], float tension = 0.0f, float overshootRatio = 0.5f, size_t stride = sizeof(Vec3f));
    // A version of SplinesFromPoints that adjusts tension per point to avoid overshoots/loops
    // caused by disparate segment sizes.

    int   SplinesFromBezier (int numPoints, const Vec3f points[], const Vec3f hullPoints[], Spline3 splines[], bool split = false); // Creates splines from the given points and Bezier hull points. If 'split' is false the splines are assumed to be continuous and numPoints - 1 splines are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.
    int   SplinesFromHermite(int numPoints, const Vec3f points[], const Vec3f tangents  [], Spline3 splines[], bool split = false); // Creates splines from the given points and tangents. If 'split' is false the splines are assumed to be continuous and numPoints - 1 splines are output. Otherwise the points are assumed to come in pairs, and numPoints / 2 splines output.

    void  MakeMonotonic(int numSplines, Spline3 splines[], bool closed = false); // closed = splines are a closed path
    // Makes splines monotonic (internal spline values always within the end points) by adjusting tangents. Useful for situations where overshoots lead to invalid values, e.g., RGB colours.

    // Queries
    Vec3f Position0(const Spline3& spline); // Start of spline
    Vec3f Position1(const Spline3& spline); // End of spline

    Vec3f Velocity0(const Spline3& spline); // Start velocity (tangent x dp/dt)
    Vec3f Velocity1(const Spline3& spline); // End velocity

    Vec3f Position    (const Spline3& spline, float t); // Returns position at 't'
    Vec3f Velocity    (const Spline3& spline, float t); // Returns velocity at 't'
    Vec3f Acceleration(const Spline3& spline, float t); // Returns acceleration at 't'
    Vec3f Jerk        (const Spline3& spline, float t); // Returns jerk at 't' -- for cubic splines this is constant over the spline

    float Curvature   (const Spline3& spline, float t); // Returns curvature at 't'. Curvature = 1 / r, where r is the radius of the osculating circle, so 0 on a straight path.
    float Torsion     (const Spline3& spline, float t); // Returns torsion at 't'. This is a measure of how the plane of curvature rotates along the path, so 0 for a flat path. Unlike curvature it is signed.

    Mat3f Frame       (const Spline3& spline, float t);               // Returns the classic Frenet frame at 't': tangent/normal/binormal. x=forward, y=left/right depending on direction of curvature, z=up/down depending on y.
    void  Frame       (const Spline3& spline, float t, Mat3f* frame); // Updates 'frame' using the Frenet frame at 't', maintaining consistency by avoiding any axis flips.
    Mat3f FrameX      (const Spline3& spline, float t);               // Returns frame with x forwards and y=left, z=up
    Mat3f FrameX      (const Spline3& spline, float t, Vec3f up);     // Returns right-handed frame with x forwards and the given up vector.

    float LengthEstimate(const Spline3& s, float* error);             // Returns estimate of length of s and optionally in 'error' the maximum error of that length.
    float Length        (const Spline3& s, float maxError = 0.01f);   // Returns length of spline accurate to the given tolerance
    float Length        (const Spline3& s, float t0, float t1, float maxError = 0.01f); // Returns length of spline segment over [t0, t1].

    Bounds3 FastBounds (const Spline3& spline); // Returns fast, convervative bounds based off the convex hull of the spline.
    Bounds3 ExactBounds(const Spline3& spline); // Returns exact bounds, taking into account extrema, requires solving a quadratic

    // Line conversion
    typedef void tEmitLinesFunc3(int numLines, Vec3f lines[][2], void* context); // Called with up to kMaxEmitLines lines at a time
    int SplinesToLinesLinear  (int numSplines, const Spline3 splines[], tEmitLinesFunc3 emitLines, void* context, float step = 0.05f); // Convert given splines into lines using even steps in 't'
    int SplinesToLinesAdaptive(int numSplines, const Spline3 splines[], tEmitLinesFunc3 emitLines, void* context, float tolerance = 0.05f); // Convert given splines into lines adaptively, with given tolerance in world units.

    typedef void tEmitLinesExFunc3(int numLines, Vec3f lines[][2], const Spline3& spline, int i, int n, float ts[][2], void* context); // Called with up to kMaxEmitLines lines at a time, plus t values and owning spline.
    int SplinesToLinesLinear  (int numSplines, const Spline3 splines[], tEmitLinesExFunc3 emitLines, void* context, float step = 0.05f);
    int SplinesToLinesAdaptive(int numSplines, const Spline3 splines[], tEmitLinesExFunc3 emitLines, void* context, float tolerance = 0.05f);

    // Subdivision
    void Split(const Spline3& spline,          Spline3* spline0, Spline3* spline1); // Splits 'spline' into two halves (at t = 0.5) and stores the results in 'subSplines'
    void Split(const Spline3& spline, float t, Spline3* spline0, Spline3* spline1); // Fills 'subSplines' with the splines corresponding to [0, t] and [t, 1]
    bool Join (const Spline3& spline0, const Spline3& spline1, Spline3* spline);    // Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    void Split(std::vector<Spline3>* splines);          // Subdivide each spline into two pieces
    void Split(std::vector<Spline3>* splines, int n);   // Subdivide each spline into 'n' pieces
    void Join (std::vector<Spline3>* splines);          // Join adjacent splines where possible

    Spline3 Trim(const Spline3& spline, float t0, float t1); // Returns subsection of 'spline' between t0 and t1

    void SubdivideForLength(std::vector<Spline3>* splines, float relativeError = 0.01f);    // Subdivide splines to be close to linear, according to relativeError.
    void SubdivideForT     (std::vector<Spline3>* splines, float error = 0.01f);            // Subdivide splines to be close to linear in t, i.e., arcLength

    void SubdivideForLengthRatio(std::vector<Vec3f>& positions, float maxRatio = 8.0f);    // Insert additional positions as necessary to ensure no segment is more than 'maxRatio' times its neighbours. Useful as preconditioner for SplinesFromPoints.

    // Nearest point
    float FindClosestPoint(Vec3f p, const Spline3& spline); // Returns t value of the closest point on s to 'p'
    float FindClosestPoint(Vec3f p, int numSplines, const Spline3 splines[], int* index); // Returns index of nearest spline, and 't' value of nearest point on that spline.

    struct SubSpline3
    {
        Spline3 mSpline;
        int     mParent;
        float   mD2;
    };

    int   FindNearbySplines(Vec3f p, int numSplines, const Spline3 splines[],       std::vector<SubSpline3>* nearbySplines, float* smallestFarOut = 0, int maxIter = 2);
    float FindClosestPoint (Vec3f p, int numSplines, const Spline3 splines[], const std::vector<SubSpline3>& nearbySplines, int* index);

    int   FindSplineIntersections(const Spline3& spline0, const Spline3& spline1, int maxResults, float results[][2], float tolerance = 0.1f);
    // Returns up to 'maxResults' intersections between the two splines, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines0, const Spline3 splines0[], int numSplines1, const Spline3 splines1[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    // Returns up to 'maxResults' intersections between the two spline lists, accurate to the given tolerance.
    int   FindSplineIntersections(int numSplines, const Spline3 splines[], int maxResults, int resultsI[][2], float resultsT[][2], float tolerance = 0.1f);
    // Returns up to 'maxResults' self-intersections in the given spline list, accurate to the given tolerance.

    // Linear point movement along spline set
    float AdvanceAgent(const Spline3& spline, float t, float dl); // Advances 'agent' at 't' on the given spline, by dl (delta length), returning t' of the new location.

    bool  AdvanceAgent(int* index, float* t, int numSplines, const Spline3 splines[], float dl); // Version of AdvanceAgent for a set of splines, assumed to be continuous. Returns false if the agent has gone off the end, in which case call ClampAgent/WrapAgent/ReverseAgent.
    void  ClampAgent  (int* index, float* t, int numSplines);    // Clamps agent position back to the nearest endpoint.
    void  WrapAgent   (int* index, float* t, int numSplines);    // Wraps the agent from the end back to the start, or vice versa.
    void  ReverseAgent(int* index, float* t);                    // Reverses the agent. (You must also negate the sign of dl yourself.)

    // Misc operations
    Spline3 Reverse(const Spline3& spline);                      // Reverses spline endpoints and tangents so that g(t) = f(1 - t).
    Spline3 Offset (const Spline3& spline, float offset);        // Offset spline, e.g., for stroking, +ve = to the right.

    void  Reverse(std::vector<Spline3>* splines);                // Reverses entire spline list
    void  Offset (std::vector<Spline3>* splines, float offset);  // Offset splines, e.g., for stroking, +ve = to the right.
}


// --- Inlines -----------------------------------------------------------------

inline int SL::NumSplinesForPoints(int numPoints)
{
    if (numPoints < 2)
        return numPoints;

    return numPoints - 1;
}

namespace
{
    inline Vec4f BezierWeights(float t)
    // Returns Bezier basis weights for 't'
    {
        float s  = 1.0f - t;

        float t2 = t * t;
        float t3 = t2 * t;

        float s2 = s * s;
        float s3 = s2 * s;

        return Vec4f(s3, 3.0f * s2 * t, 3.0f * s * t2, t3);
    }

    inline Vec4f BezierWeightsD1(float t)
    {
        float t2 = t * t;

        return Vec4f
        (
            - 3.0f +  6.0f * t - 3.0f * t2,
              3.0f - 12.0f * t + 9.0f * t2,
                      6.0f * t - 9.0f * t2,
                                 3.0f * t2
        );
    }

    inline Vec4f BezierWeightsD2(float t)
    {
        return Vec4f
        (
             6.0f -  6.0f * t,
           -12.0f + 18.0f * t,
             6.0f - 18.0f * t,
                     6.0f * t
        );
    }

    const Vec4f kBezierWeightsD3(-6.0f, 18.0f, -18.0f, 6.0f);
}


// 1D

inline float SL::Position0(const Spline1& spline)
{
    return spline.x;
}
inline float SL::Position1(const Spline1& spline)
{
    return spline.w;
}

inline float SL::Velocity0(const Spline1& spline)
{
    return 3.0f * (spline.y - spline.x);
}
inline float SL::Velocity1(const Spline1& spline)
{
    return 3.0f * (spline.w - spline.z);
}

inline float SL::Position(const Spline1& spline, float t)
{
    return dot(spline, BezierWeights(t));
}

inline float SL::Velocity(const Spline1& spline, float t)
{
    return dot(spline, BezierWeightsD1(t));
}

inline float SL::Acceleration(const Spline1& spline, float t)
{
    return dot(spline, BezierWeightsD2(t));
}

inline float SL::Jerk(const Spline1& spline, float)
{
    return dot(spline, kBezierWeightsD3);
}

inline Vec4f SL::CubicCoeffs(const Spline1& b)
{
    return Vec4f
    (
                b.x                                ,
        -3.0f * b.x + 3.0f * b.y                   ,
         3.0f * b.x - 6.0f * b.y + 3.0f * b.z      ,
               -b.x + 3.0f * b.y - 3.0f * b.z + b.w
    );
}

// 2D

inline Vec2f SL::Position0(const Spline2& spline)
{
    return Vec2f(spline.xb.x, spline.yb.x);
}
inline Vec2f SL::Position1(const Spline2& spline)
{
    return Vec2f(spline.xb.w, spline.yb.w);
}

inline Vec2f SL::Velocity0(const Spline2& spline)
{
    return 3.0f * Vec2f
    (
        spline.xb.y - spline.xb.x,
        spline.yb.y - spline.yb.x
    );
}
inline Vec2f SL::Velocity1(const Spline2& spline)
{
    return 3.0f * Vec2f
    (
        spline.xb.w - spline.xb.z,
        spline.yb.w - spline.yb.z
    );
}

inline void SL::ClampAgent(int* index, float* t, int numSplines)
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

inline void SL::WrapAgent(int* indexInOut, float* tInOut, int numSplines)
{
    int& index = *indexInOut;
    float& t = *tInOut;

    SL_ASSERT(!IsNAN(t));
    SL_ASSERT(index == 0 || index == numSplines - 1);

    t -= floorf(t);
    index ^= numSplines - 1;
}

inline void SL::ReverseAgent(int* , float* t)
{
    *t = ceilf(*t) - *t;
}

// 3D

inline Vec3f SL::Position0(const Spline3& spline)
{
    return Vec3f(spline.xb.x, spline.yb.x, spline.zb.x);
}
inline Vec3f SL::Position1(const Spline3& spline)
{
    return Vec3f(spline.xb.w, spline.yb.w, spline.zb.w);
}

inline Vec3f SL::Velocity0(const Spline3& spline)
{
    return 3.0f * Vec3f
    (
        spline.xb.y - spline.xb.x,
        spline.yb.y - spline.yb.x,
        spline.zb.y - spline.zb.x
    );
}
inline Vec3f SL::Velocity1(const Spline3& spline)
{
    return 3.0f * Vec3f
    (
        spline.xb.w - spline.xb.z,
        spline.yb.w - spline.yb.z,
        spline.zb.w - spline.zb.z
    );
}

#endif
