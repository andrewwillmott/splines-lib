//
//  RotSplines.hpp
//
//  Cubic spline utilities for rotations
//
//  Andrew Willmott
//

#ifndef SL_ROT_SPLINES_H
#define SL_ROT_SPLINES_H

#include "Splines.hpp"

namespace SL
{
    // 2D
    typedef Spline1 RotSpline2;

    RotSpline2 BezierRotSpline (float r0, float r1, float r2, float r3);    // Rotation angles in radians
    RotSpline2 HermiteRotSpline(float r0, float r1, float wv0, float wv1);  // Start/end rotations + rot velocity in/out
    RotSpline2 LineRotSpline   (float r0, float r1);                        // Simple linear interpolation

    float Rotation0(const RotSpline2& rs);  // Start rotation angle
    float Rotation1(const RotSpline2& rs);  // End rotation angle

    float RotVelocity0(const RotSpline2& rs);  // Start rotational velocity (signed rad/s)
    float RotVelocity1(const RotSpline2& rs);  // End rotational velocity (signed rad/s)

    float Rotation       (const RotSpline2& rs, float t);  // Interpolated rotation angle. Use with RRot2f() to get corresponding matrix.
    float RotVelocity    (const RotSpline2& rs, float t);  // Corresponding rotational velocity (rad/s).
    float RotAcceleration(const RotSpline2& rs, float t);  // Corresponding rotational acceleration (rad/s/s).

    float AngSpeed       (const RotSpline2& rs, float t);  // Positive angular rotation (mostly for consistency with 3D api below)
    float AngAcceleration(const RotSpline2& rs, float t);  // Positive angular acceleration

    int RotSplinesFromAngles(int numAngles, const float a[], RotSpline2 splines[], float tension = 0.0f, size_t stride = sizeof(float));
    int RotSplinesFromDirs  (int numDirs  , const Vec2f d[], RotSpline2 splines[], float tension = 0.0f, size_t stride = sizeof(Vec2f));
    // Creates a series of RotSplines that smoothly interpolate the given set of unit quaternions, parameterised by 'tension'.

    // 3D
    struct RotSpline3
    {
        Quatf q0;
        Vec3f w1;  // Log(q0 -> q1)
        Vec3f w2;  // Log(q1 -> q2)
        Vec3f w3;  // Log(q2 -> q3)
    };

    RotSpline3 BezierRotSpline (Quatf q0, Quatf q1, Quatf q2, Quatf q3);    // q1/q2 are the usual control points
    RotSpline3 HermiteRotSpline(Quatf q0, Quatf q1, Vec3f w0, Vec3f w1);    // w0/w1 here are rot velocity, e.g., from QuatDiff3(q-1, q+1)
    RotSpline3 LineRotSpline   (Quatf q0, Quatf q1);  // SLerp equivalent

    Quatf Rotation0(const RotSpline3& rs);  // Start rotation
    Quatf Rotation1(const RotSpline3& rs);  // End rotation

    Vec3f RotVelocity0(const RotSpline3& rs);  // Start rotational velocity
    Vec3f RotVelocity1(const RotSpline3& rs);  // End rotational velocity

    // Smooth rotation interpolation using 'cumulative' approach of KKS, with well-defined cheap(ish) derivatives. C2.
    Quatf Rotation       (const RotSpline3& rs, float t);  // Interpolated unit quaternion rotation. Use with RRot3f() to get corresponding matrix.
    Vec3f RotVelocity    (const RotSpline3& rs, float t);  // Derivative of Rotation() in tangent space (rotation axis x dtheta/dt).
    Quatf RotAcceleration(const RotSpline3& rs, float t);  // Second derivative of Rotation() in tangent space (axis x d2theta/dt2 x scale).

    // These scalar angle variants match their 2D variants. Useful for comparison or you're interested in a measure of magnitude.
    float AngSpeed       (const RotSpline3& rs, float t);  // Angular speed: 2 x len(rv), with the 2 coming from the double cover.
    float AngAcceleration(const RotSpline3& rs, float t);  // Angular acceleration: 2 x len(xyz(ra)). (Picking theta out of the log form of the derivative.)

    Quatf Rotation(Quatf q0, Quatf q1, Quatf q2, Quatf q3, float t);     // Direct rotation interpolation using classic four Bezier-form quats
    Quatf Rotation(Quatf q0, Vec3f ld1, Vec3f ld2, Vec3f ld3, float t);  // More efficient rotation interpolation using initial quat and three 1/3 log-form deltas

    // Alternate 'Shoemake' form of interpolation, using SLerp + De Casteljau algorithm. Not as amenable to derivatives/analysis. C1.
    Quatf RotationShoemake(const RotSpline3& rs, float t);
    Quatf RotationShoemake(Quatf q0, Quatf q1, Quatf q2, Quatf q3, float t);

    // Alternate 'Spherical Quadrangle' interpolation: simpler/faster version of Shoemake that replaces three linear slerps with one quadratic. C1 but notably jerkier at end points.
    Quatf RotationSQuad(const RotSpline3& rs, float t);
    Quatf RotationSQuad(Quatf q0, Quatf q1, Quatf q2, Quatf q3, float t);  // Note: technically need to set q1/q2 with 1/2 factor rather than 1/3, but not found much difference.

    void Split(const RotSpline3& rs, float t, RotSpline3* rs0, RotSpline3* rs1);    // Splits 'rs' into two halves (at t = 0.5) and stores the results in rs0/1
    bool Join (const RotSpline3& rs0, const RotSpline3& rs1, RotSpline3* rs);       // Joins two splines that were formerly Split(). Assumes t=0.5, returns false if the source splines don't match up.

    int RotSplinesFromQuats(int numQuats, const Quatf qi[], RotSpline3 splines[], float tension = 0.0f, size_t stride = sizeof(Quatf));
    // Creates a series of RotSplines that smoothly interpolate the given set of unit quaternions, parameterised by 'tension'.
}


// Inlines

// would be nice if they ever bothered to add a robust function alias
inline SL::RotSpline2 SL::BezierRotSpline (float w0, float w1, float w2, float w3)
{
    return BezierSpline(w0, w1, w2, w3);
}

inline SL::RotSpline2 SL::HermiteRotSpline(float w0, float w1, float wv0, float wv1)
{
    return HermiteSpline(w0, w1, wv0, wv1);
}

inline SL::RotSpline2 SL::LineRotSpline(float w0, float w1)
{
    return LineSpline(w0, w1);
}

inline float SL::Rotation0(const RotSpline2& rs)
{
    return Position0(rs);
}

inline float SL::Rotation1(const RotSpline2& rs)
{
    return Position1(rs);
}

inline float SL::RotVelocity0(const RotSpline2& rs)
{
    return Velocity0(rs);
}

inline float SL::RotVelocity1(const RotSpline2& rs)
{
    return Velocity1(rs);
}

inline float SL::Rotation(const RotSpline2& rs, float t)
{
    return Position(rs, t);
}

inline float SL::RotVelocity(const RotSpline2& rs, float t)
{
    return Velocity(rs, t);
}

inline float SL::RotAcceleration(const RotSpline2& rs, float t)
{
    return Acceleration(rs, t);
}

inline float SL::AngSpeed(const RotSpline2& rs, float t)
{
    return fabsf(Velocity(rs, t));
}

inline float SL::AngAcceleration(const RotSpline2& rs, float t)
{
    return fabsf(Acceleration(rs, t));
}

inline float SL::AngSpeed(const RotSpline3& rs, float t)
{
    return 2 * len(RotVelocity(rs, t));
}

inline float SL::AngAcceleration(const RotSpline3& rs, float t)
{
    return 2 * len(xyz(RotAcceleration(rs, t)));
}

#endif
