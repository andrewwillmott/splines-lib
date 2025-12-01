//
// RotSplines.cpp
//
// Cubic spline utilities for rotations
//
// Andrew Willmott
//

#include "RotSplines.hpp"

using namespace SL;

namespace
{
    inline float ConstrainAngle(float ra, float rb)
    {
        float d = (rb - ra) * (1.0f / vlf_twoPi);
        float id = floorf(d + 0.5f);
        return ra + (d - id) * vlf_twoPi;
    }

    inline RotSpline2 RotSplineFromAngles(const float a[], int i0, int i1, int i2, int i3, float s)
    {
        float pb1 = a[i1] + s * (a[i2] - a[i0]);
        float pb2 = a[i2] - s * (a[i3] - a[i1]);

        return BezierSpline(a[i1], pb1, pb2, a[i2]);
    }
}

int SL::RotSplinesFromAngles(int numAngles, const float ai[], RotSpline2 splines[], float tension, size_t stride)
{
    SL_ASSERT(numAngles >= 0);

    float s = (1.0f - tension) * (1.0f / 6.0f);  // 1/2 for averaging * 1/3 for v scale

    float local[4];
    int ln = vl_min(numAngles, 4);

    for (int i = 0; i < ln; i++)
    {
        local[i] = *ai;
        (char*&) ai += stride;
    }

    for (int i = 1; i < ln; i++)
        local[i] = ConstrainAngle(local[i - 1], local[i]);

    switch (numAngles)
    {
    case 0: return 0;
    case 1: *splines = RotSplineFromAngles(local, 0, 0, 0, 0, s); return 1;
    case 2: *splines = RotSplineFromAngles(local, 0, 0, 1, 1, s); return 1;
    }

    *splines++ = RotSplineFromAngles(local, 0, 0, 1, 2, s);
    int base = 0;
    int i0 = 0, i1 = 1, i2 = 2, i3 = 3;

    while (true)
    {
        *splines++ = RotSplineFromAngles(local, i0, i1, i2, i3, s);

        i0 = i1; i1 = i2; i2 = i3; i3 = (base++ & 3);

        if (base == numAngles - 3)
            break;

        local[i3] = ConstrainAngle(local[i2], *ai);
        (char*&) ai += stride;
    }

    *splines++ = RotSplineFromAngles(local, i0, i1, i2, i2, s);

    return numAngles - 1;
}

int SL::RotSplinesFromDirs(int numAngles, const Vec2f d[], RotSpline2 splines[], float tension, size_t stride)
{
    SL_ASSERT(numAngles >= 0);

    float s = (1.0f - tension) * (1.0f / 6.0f);  // 1/2 for averaging * 1/3 for v scale

    float local[4];
    int ln = vl_min(numAngles, 4);

    for (int i = 0; i < ln; i++)
    {
        local[i] = atan2f(d->y, d->x);
        (char*&) d += stride;
    }

    for (int i = 1; i < ln; i++)
        local[i] = ConstrainAngle(local[i - 1], local[i]);

    switch (numAngles)
    {
    case 0: return 0;
    case 1: *splines = RotSplineFromAngles(local, 0, 0, 0, 0, s); return 1;
    case 2: *splines = RotSplineFromAngles(local, 0, 0, 1, 1, s); return 1;
    }

    *splines++ = RotSplineFromAngles(local, 0, 0, 1, 2, s);
    int base = 0;
    int i0 = 0, i1 = 1, i2 = 2, i3 = 3;

    while (true)
    {
        *splines++ = RotSplineFromAngles(local, i0, i1, i2, i3, s);

        i0 = i1; i1 = i2; i2 = i3; i3 = (base++ & 3);

        if (base == numAngles - 3)
            break;

        local[i3] = ConstrainAngle(local[i2], atan2f(d->y, d->x));
        (char*&) d += stride;
    }

    *splines++ = RotSplineFromAngles(local, i0, i1, i2, i2, s);

    return numAngles - 1;
}

//
// Note: if adapting this code for use with an alternative quaternion library,
// be aware that to make the code more readable, the quaternion multiplication
// routines here work left-to-right. I.e., QuatMult(a, b) = apply rotation a
// and then rotation b. QuatApply() is the standard transformation of a 3d
// point by a vector, which has no such ordering issues.
//

namespace
{
    // b   = 1,  1 - (1 - t)^3,  3t^2 - 2t^3,    t^3
    // b'  = 0,     3(1 - t)^2,  6t   - 6t^2,  3 t^2
    // b'' = 0,    -6(1 - t)  ,  6    - 12t ,  6 t

    inline Vec3f BezierRotWeights(float t)
    {
        float t2 = sqr(t);
        float t3 = t2 * t;
        float it = 1 - t;
        float it3 = it * it * it;

        return { 1 - it3, 3 * t2 - 2 * t3, t3 };
    }

    inline Vec3f BezierRotWeightsD1(float t)
    {
        float t2 = sqr(t);
        float it = 1.0f - t;
        float it2 = sqr(it);

        return { 3 * it2, 6 * t * it, 3 * t2 };
    }

    inline Vec3f BezierRotWeightsD2(float t)
    {
        return { -6 + 6 * t, 6 - 12 * t, 6 * t };
    }
}

RotSpline3 SL::BezierRotSpline(Quatf q0, Quatf q1, Quatf q2, Quatf q3)
{
    return RotSpline3
    {
        q0,
        QuatDiff3(q0, q1),
        QuatDiff3(q1, q2),
        QuatDiff3(q2, q3)
    };
}

#if 0
RotSpline3 SL::HermiteRotSplineRef(Quatf q0, Quatf q1, Vec3f w0, Vec3f w1)
{
    w0 *= float(1.0 / 3.0);
    w1 *= float(1.0 / 3.0);

    // Easier to understand Bezier form, following
    //   qm0 = q0 + w0 / 3
    //   qm1 = q0 - w1 / 3
    Quatf qm0 = QuatMult(q0, ExpUnit3(+w0));
    Quatf qm1 = QuatMult(q1, ExpUnit3(-w1)); // q1 exp(w1)-1 = q1 exp(-w1)
    return BezierRotSpline(q0, qm0, qm1, q1);
}
#endif

RotSpline3 SL::HermiteRotSpline(Quatf q0, Quatf q1, Vec3f w0, Vec3f w1)
{
    w0 *= float(1.0 / 3.0);
    w1 *= float(1.0 / 3.0);

    // More efficient: find middle delta 'wm' directly and re-use w0/w1
    // q0 e^w0 e^wm e^w1 = q1 -> e^wm = e^-w0 q0-1 q1 e^-w1
    Quatf qm;
    qm =              ExpUnit3(-w0);
    qm = QuatMult(qm, QuatMult(QuatConj(q0), q1));
    qm = QuatMult(qm, ExpUnit3(-w1));

    Vec3f wm = LogUnit3(qm);

    return RotSpline3{ q0, w0, wm, w1 };
}

RotSpline3 SL::LineRotSpline(Quatf q0, Quatf q1)
{
    RotSpline3 rs;
    rs.q0 = q0;

    Vec3f w = QuatDiff3(q0, q1) * (1.0f / 3.0f);
    rs.w1 = w;
    rs.w2 = w;
    rs.w3 = w;

    return rs;
}

Quatf SL::Rotation0(const RotSpline3& rs)
{
    return rs.q0;
}

Quatf SL::Rotation1(const RotSpline3& rs)
{
    Quatf r = rs.q0;
    r = QuatMult(r, ExpUnit3(rs.w1));
    r = QuatMult(r, ExpUnit3(rs.w2));
    r = QuatMult(r, ExpUnit3(rs.w3));
    return r;
}

Vec3f SL::RotVelocity0(const RotSpline3& rs)
{
    return rs.w1 * 3;
}

Vec3f SL::RotVelocity1(const RotSpline3& rs)
{
    return rs.w3 * 3;
}

Quatf SL::Rotation(const RotSpline3& rs, float t)
{
    // Note: in contrast to the geometric construction in Shoemake's approaches,
    // this routine uses Kim/Kim/Shin's general algebraic construction scheme
    // for unit quaternion splines:
    //
    //   https://www.semanticscholar.org/paper/A-general-construction-scheme-for-unit-quaternion-Kim-Kim/95336585fc493e66ed79132d2fefd4d157f8cc49
    //
    // The basic idea is to use the exponential map to do the basis function
    // interpolation in axis/angle space. If a unit quaternion q is written as
    // [n sin(w), cos(w)], then log(q) = [n w, 0], and the exponent of that gets
    // you back to the original quaternion as you might expect.
    //
    // In 2D, the equivalent of this, as implemented above, is to operate
    // directly on angles, and then at the end convert the result to an
    // orientation vector via sin and cos (the 2D quaternion equivalent), and
    // hence to a matrix. You may be familiar with the corresponding exp/log
    // relationship: exp(i w) = i sin(w), cos(w).
    //
    // In 3D, because rotations don't commute, the addition part of this
    // approach no longer works, because exp(a) exp(b) != exp(a + b). (There is
    // a formula for exp(a) exp(b) in terms of AB - BA differences, but it does
    // not converge well in the case of quaternions.) We can still apply the 't'
    // part of 'a + t(b - a)' (say) in log(q) space, but have to transform the
    // result back to quaternion space to concatenate the parts. Hence the
    // classic log form of SLerp isn't
    //
    //   q = exp(log(a) + t (log(b) - log(a)))
    //
    // but rather
    //
    //   q = a exp(t log(b-1 a))
    //
    // For cubic splines, we take the same approach, but now we must concatenate
    // four of these terms (one for each basis weight/control point), and,
    // because we are not free to re-arrange terms due to commutivity, the
    // standard Bezier weights must be recast in cumulative form.
    //
    // Conceptually this is a little more complex than the geometric
    // formulations that utilise nested slerps, but with the right choice of
    // tangents (namely in/out rotation velocities) the results can be C2
    // continuous, and, the log/exp framework makes derivatives straightforward.

    // Right, let's do this thing. The code itself is pretty straight-forward.
    Vec3f b = BezierRotWeights(t);  // cumulative weights rather than individual BezierWeights(). First implicit term is thus always 1.

    Quatf r = rs.q0;
    r = QuatMult(r, ExpUnit3(b.x * rs.w1));  // Basically three quat multiplies + 3 sin/cos + 3v*. No acos/atan.
    r = QuatMult(r, ExpUnit3(b.y * rs.w2));  // Shoemake original approach is 6 slerps. SQuad is 3 slerps.
    r = QuatMult(r, ExpUnit3(b.z * rs.w3));

    return r;
}

// #define USE_REFERENCE

namespace
{
    template <typename... T> inline Quatf QuatMult(const Quatf& q0, const Quatf& q1, T... args)
    {
        return QuatMult(QuatMult(q0, q1), args...);
    }

#ifdef USE_REFERENCE
    Vec3f RotVelocityRef(const RotSpline3& rs, float t)
    {
        // Everything written out for clarity. We're implementing:
        // q  = q0 exp(w1 b1) exp(w2 b2) exp(w3 b3)
        // q' =
        //   q0 exp(w1 b1) (w1 b'1) exp(w2 b2) exp(w3 b3)
        // + q0 exp(w1 b1) exp(w2 b2) (w2 b'2) exp(w3 b3)
        // + q0 exp(w1 b1) exp(w2 b2) exp(w3 b3) (w3 b'3)
        //
        // The above looks more complex than the linear form that can be used
        // in the various other Velocity routines, but it's just the chain
        // rule applied to the multiplied terms, vs the additions elsewhere.

        Vec3f b0 = BezierRotWeights  (t);
        Vec3f b1 = BezierRotWeightsD1(t);

        Quatf qd1(ExpUnit3(rs.w1 * b0.x));
        Quatf qd2(ExpUnit3(rs.w2 * b0.y));
        Quatf qd3(ExpUnit3(rs.w3 * b0.z));

        Quatf qv1; xyz(qv1) = rs.w1 * b1.x; qv1.w = 0;
        Quatf qv2; xyz(qv2) = rs.w2 * b1.y; qv2.w = 0;
        Quatf qv3; xyz(qv3) = rs.w3 * b1.z; qv3.w = 0;

        Quatf qv;
        qv  = QuatMult(rs.q0, qd1, qv1, qd2, qd3);
        qv += QuatMult(rs.q0, qd1, qd2, qv2, qd3);
        qv += QuatMult(rs.q0, qd1, qd2, qd3, qv3);

        // At this point we have q'(t) = q(t) w'(t), analogous to (e^iw)' = w' e^iw
        // so we need to invert out the q(t) part to get w'(t)
        Quatf q = QuatMult(rs.q0, qd1, qd2, qd3);
        Quatf w = QuatMult(QuatConj(q), qv);

        return xyz(w);
    }
#endif
}

Vec3f SL::RotVelocity(const RotSpline3& rs, float t)
{
    // This is a much-simplified version of RotVelocityRef with terms cancelled out.
    // 1 QuatMult, 2 QuatApply, 2 ExpUnit3=sincos, vs Rotation = 3 QM, 3 EU
    Vec3f b = BezierRotWeights  (t);
    Vec3f v = BezierRotWeightsD1(t);

    Quatf qd2  = ExpUnit3(rs.w2 * b.y);
    Quatf qd3  = ExpUnit3(rs.w3 * b.z);

    Quatf qd23  = QuatMult(qd2, qd3);

    Vec3f wv1(rs.w1 * v.x);
    Vec3f wv2(rs.w2 * v.y);
    Vec3f wv3(rs.w3 * v.z);

    Vec3f w;
    w  = QuatApply(wv1, qd23);
    w += QuatApply(wv2, qd3 );
    w +=           wv3;

#ifdef USE_REFERENCE
    SL_ASSERT(len(w - RotVelocityRef(rs, t)) < 1e-5f);
#endif

    return w;
}

namespace
{
    // Quat/Vector Quat multiply
    inline Quatf QuatMult(const Quatf& a, const Vec3f& b)
    {
        Quatf q;
        xyz(q) = cross(xyz(a), b) + a.w * b;
        q.w    =  -dot(xyz(a), b);
        return q;
    }
}

Quatf SL::RotAcceleration(const RotSpline3& rs, float t)
{
    /*
        A01        = exp(w1 b1)
        A02        = exp(w2 b2)
        A03        = exp(w3 b3)

        A11 = A01' = A01 (w1 b'1)
        A12 = A02' = A02 (w2 b'2)
        A13 = A03' = A03 (w3 b'3)

        A21 = A11' = A01 (w1 b''1) + A11 (w1 b'1)
        A22 = A12' = A02 (w2 b''1) + A12 (w2 b'1)
        A23 = A13' = A03 (w3 b''1) + A13 (w3 b'1)

        q' =      (A11 A02 A03 + A01 A12 A03 + A01 A02 A13)

        q'' =     (A21 A02 A03 + A01 A22 A03 + A01 A02 A23)
            + 2 * (A11 A12 A03 + A11 A02 A13 + A01 A12 A13)
    */

    Vec3f b0 = BezierRotWeights  (t);
    Vec3f b1 = BezierRotWeightsD1(t);
    Vec3f b2 = BezierRotWeightsD2(t);

    Quatf a01(ExpUnit3(rs.w1 * b0.x));
    Quatf a02(ExpUnit3(rs.w2 * b0.y));
    Quatf a03(ExpUnit3(rs.w3 * b0.z));

    Vec3f wv1(rs.w1 * b1.x);
    Vec3f wv2(rs.w2 * b1.y);
    Vec3f wv3(rs.w3 * b1.z);

    Vec3f wa1(rs.w1 * b2.x);
    Vec3f wa2(rs.w2 * b2.y);
    Vec3f wa3(rs.w3 * b2.z);

    Quatf a11 = QuatMult(a01, wv1);
    Quatf a12 = QuatMult(a02, wv2);
    Quatf a13 = QuatMult(a03, wv3);

    Quatf a21 = QuatMult(a11, wv1) + QuatMult(a01, wa1); // a01 (wv1 wv1 + wa1)
    Quatf a22 = QuatMult(a12, wv2) + QuatMult(a02, wa2);
    Quatf a23 = QuatMult(a13, wv3) + QuatMult(a03, wa3);

    Quatf wd = QuatMult(a01, a02, a03);

    Quatf qa1;
    qa1  = QuatMult(a21, a02, a03);
    qa1 += QuatMult(a01, a22, a03);
    qa1 += QuatMult(a01, a02, a23);

    Quatf qa2;
    qa2  = QuatMult(a11, a12, a03);
    qa2 += QuatMult(a11, a02, a13);
    qa2 += QuatMult(a01, a12, a13);

    Quatf qa = qa1 + 2 * qa2;
    Quatf wa = QuatMult(QuatConj(wd), qa);

    return wa;
}

Quatf SL::Rotation(Quatf q0, Vec3f w1, Vec3f w2, Vec3f w3, float t)
{
    Vec3f b = BezierRotWeights(t);

    Quatf r = q0;
    r = QuatMult(r, ExpUnit3(b.x * w1));
    r = QuatMult(r, ExpUnit3(b.y * w2));
    r = QuatMult(r, ExpUnit3(b.z * w3));
    return r;
}

Quatf SL::Rotation(Quatf q0, Quatf q1, Quatf q2, Quatf q3, float t)
{
    Vec3f w1 = QuatDiff3(q0, q1);
    Vec3f w2 = QuatDiff3(q1, q2);
    Vec3f w3 = QuatDiff3(q2, q3);

    return Rotation(q0, w1, w2, w3, t);
}

Quatf SL::RotationShoemake(const RotSpline3& rs, float t)
{
    Quatf q0 = rs.q0;

    // We're saving the acos here, but only in the first round.
    Quatf q01 = QuatMult(q0, ExpUnit3(rs.w1 * t));
    Quatf q1  = QuatMult(q0, ExpUnit3(rs.w1));
    Quatf q12 = QuatMult(q1, ExpUnit3(rs.w2 * t));
    Quatf q2  = QuatMult(q1, ExpUnit3(rs.w2));
    Quatf q23 = QuatMult(q2, ExpUnit3(rs.w3 * t));

    Quatf q012 = SLerp(q01, q12, t);
    Quatf q123 = SLerp(q12, q23, t);

    Quatf q0123 = SLerp(q012, q123, t);

    return q0123;
}

Quatf SL::RotationShoemake(Quatf q0, Quatf q1, Quatf q2, Quatf q3, float t)
{
    // Basically, just De Casteljau's algorithm using slerps rather than linear interpolation.
    Quatf q01 = SLerp(q0, q1, t);
    Quatf q12 = SLerp(q1, q2, t);
    Quatf q23 = SLerp(q2, q3, t);

    Quatf q012 = SLerp(q01, q12, t);
    Quatf q123 = SLerp(q12, q23, t);

    Quatf q0123 = SLerp(q012, q123, t);

    return q0123;
}

Quatf SL::RotationSQuad(const RotSpline3& rs, float t)
{
    // Cheaper version where you slerp the end points and the interior points,
    // then quadratically blend between the end point result and the interior
    // result, such that you get the interior point in the middle and exterior
    // points at either end. *waves hands*.
    //
    // RotSpline3 is not the best format here, q0, w123, q1, w2 would be better,
    // but this is just for cross-checking. Use the standalone variant if this
    // is being used in anger.
    Quatf q0 = rs.q0;

    Quatf q1  = QuatMult(q0, ExpUnit3(rs.w1));
    Quatf q12 = QuatMult(q1, ExpUnit3(rs.w2 * t));

    Quatf q3  = QuatMult(q1, QuatMult(ExpUnit3(rs.w2), ExpUnit3(rs.w3)));  // :|
    Quatf q03 = SLerp(q0, q3, t);

    Quatf q0123 = SLerp(q03, q12, 2 * t * (1 - t));

    return q0123;
}

Quatf SL::RotationSQuad(Quatf q0, Quatf q1, Quatf q2, Quatf q3, float t)
{
    Quatf q03 = SLerp(q0, q3, t);
    Quatf q12 = SLerp(q1, q2, t);

    Quatf q0123 = SLerp(q03, q12, 2 * t * (1 - t));

    return q0123;
}

void SL::Split(const RotSpline3& rs, float t, RotSpline3* spline0, RotSpline3* spline1)
{
    Quatf qm = Rotation(rs, t);
    Vec3f wm = RotVelocity(rs, t);

    Quatf q0 = Rotation0(rs);
    Quatf q1 = Rotation1(rs);

    Vec3f w0 = RotVelocity0(rs);
    Vec3f w1 = RotVelocity1(rs);

    *spline0 = HermiteRotSpline(q0, qm, w0, wm);
    *spline1 = HermiteRotSpline(qm, q1, wm, w1);
}

bool SL::Join(const RotSpline3& rs0, const RotSpline3& rs1, RotSpline3* rs)
{
    if (rs0.w3 != rs0.w1) // early out
        return false;

    *rs = HermiteRotSpline(Rotation0(rs0), Rotation1(rs1), RotVelocity0(rs0), RotVelocity1(rs1));
    return true;
}

namespace
{
    inline Quatf ConstrainQuat(const Quatf& q1, const Quatf& q2)
    {
        return (dot(q1, q2) < 0) ? -q2 : q2;
    }

    inline RotSpline3 RotSplineFromQuats(const Quatf q[], int i0, int i1, int i2, int i3, float s)
    {
        Vec3f wb1 = s * QuatDiff3(q[i0], q[i2]);
        Vec3f wb2 = s * QuatDiff3(q[i1], q[i3]);

        return HermiteRotSpline(q[i1], q[i2], wb1, wb2);
    }
}

int SL::RotSplinesFromQuats(int numQuats, const Quatf q[], RotSpline3 splines[], float tension, size_t stride)
{
    SL_ASSERT(numQuats >= 0);

    float s = (1.0f - tension) * (1.0f / 6.0f);  // 1/2 for averaging * 1/3 for v scale

    Quatf local[4];
    int ln = vl_min(numQuats, 4);

    for (int i = 0; i < ln; i++)
    {
        local[i] = *q;
        (char*&) q += stride;
    }

    for (int i = 1; i < ln; i++)
        local[i] = ConstrainQuat(local[i - 1], local[i]);

    switch (numQuats)
    {
    case 0: return 0;
    case 1: *splines = RotSplineFromQuats(local, 0, 0, 0, 0, s); return 1;
    case 2: *splines = RotSplineFromQuats(local, 0, 0, 1, 1, s); return 1;
    }

    *splines++ = RotSplineFromQuats(local, 0, 0, 1, 2, s);
    int base = 0;
    int i0 = 0, i1 = 1, i2 = 2, i3 = 3;

    while (true)
    {
        *splines++ = RotSplineFromQuats(local, i0, i1, i2, i3, s);

        i0 = i1; i1 = i2; i2 = i3; i3 = (base++ & 3);

        if (base == numQuats - 3)
            break;

        local[i3] = ConstrainQuat(local[i2], *q);
        (char*&) q += stride;
    }

    *splines++ = RotSplineFromQuats(local, i0, i1, i2, i2, s);

    return numQuats - 1;
}
