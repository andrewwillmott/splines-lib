//
// Splines.cpp
//
// Cubic spline utilities for 1D/2D/3D
//
// Andrew Willmott
//

#include "Splines.hpp"

#include <float.h>

using namespace SL;

using std::vector;

#define SL_ASSERT_INDEX(M_I, M_N)    \
    SL_ASSERT((unsigned int)(M_I) < (unsigned int)(M_N))

namespace
{
    inline float Clamp(float x, float minX, float maxX)
    {
        if (x < minX)
            return minX;
        if (x > maxX)
            return maxX;
        return x;
    }

    inline float InvSqrtFast(float x)
    {
        float xhalf = 0.5f * x;
        int32_t i = (int32_t&) x;

        i = 0x5f375a86 - (i >> 1);
        x = (float&) i;
        x = x * (1.5f - xhalf * x * x);

        return x;
    }

    inline bool Larger(const Bounds2& bb, float t) { Vec2f d = bb.mMax - bb.mMin; return d.x > t || d.y > t; }
    inline bool Larger(const Bounds3& bb, float t) { Vec3f d = bb.mMax - bb.mMin; return d.x > t || d.y > t || d.z > t; }

    inline bool Intersects(const Bounds2& a, const Bounds2& b)
    {
        return a.mMax.x >= b.mMin.x && a.mMin.x <= b.mMax.x
            && a.mMax.y >= b.mMin.y && a.mMin.y <= b.mMax.y;
    }
    inline bool Intersects(const Bounds3& a, const Bounds3& b)
    {
        return a.mMax.x >= b.mMin.x && a.mMin.x <= b.mMax.x
            && a.mMax.y >= b.mMin.y && a.mMin.y <= b.mMax.y
            && a.mMax.z >= b.mMin.z && a.mMin.z <= b.mMax.z;
    }
    inline void Add(Bounds2& a, const Bounds2& b)
    {
        if (a.mMin.x > b.mMin.x) a.mMin.x = b.mMin.x; else if (a.mMax.x < b.mMax.x) a.mMax.x = b.mMax.x;
        if (a.mMin.y > b.mMin.y) a.mMin.y = b.mMin.y; else if (a.mMax.y < b.mMax.y) a.mMax.y = b.mMax.y;
    }
    inline void Add(Bounds3& a, const Bounds3& b)
    {
        if (a.mMin.x > b.mMin.x) a.mMin.x = b.mMin.x; else if (a.mMax.x < b.mMax.x) a.mMax.x = b.mMax.x;
        if (a.mMin.y > b.mMin.y) a.mMin.y = b.mMin.y; else if (a.mMax.y < b.mMax.y) a.mMax.y = b.mMax.y;
        if (a.mMin.z > b.mMin.z) a.mMin.z = b.mMin.z; else if (a.mMax.z < b.mMax.z) a.mMax.z = b.mMax.z;
    }

    template<class T> inline int size_i(const T& container) { return int(container.size()); }
}

namespace
{
    // Utilities
    inline Vec2f ArcError2(Vec4f s)
    // Returns squared displacement from linear (b0_b3) for hull points b1/b2
    {
        float w = s.w - s.x;

        float ty = s.x + w * 1.0f / 3.0f - s.y;
        float tz = s.x + w * 2.0f / 3.0f - s.z;
        float d2 = 1.0f / (sqr(w) + 1.0f);

        return Vec2f(sqr(ty) * d2, sqr(tz) * d2);
    }

    template<class T, class S> inline S SplineFromPointsT(const char* p8, size_t stride, int i0, int i1, int i2, int i3, float tension)
    {
        T p0 = *(T*) (p8 + i0 * stride);
        T p1 = *(T*) (p8 + i1 * stride);
        T p2 = *(T*) (p8 + i2 * stride);
        T p3 = *(T*) (p8 + i3 * stride);

        float s = (1.0f - tension) * (1.0f / 6.0f);

        T pb1 = p1 + s * (p2 - p0);
        T pb2 = p2 - s * (p3 - p1);

        return BezierSpline(p1, pb1, pb2, p2);
    }

    template<class T, class S> inline S SplineFromPointsDynamicT(const char* p8, size_t stride, int i0, int i1, int i2, int i3, float tension, float ratio = 0.666f)
    {
        // The standard algorithm above has an issue when a small segment is followed by a very large one.
        // ||pb1|| can wind up being much larger than ||p0_p1||, which leads to the spline rep looping
        // back on itself. Basically, the shared tangent's length is much larger than the first segment's size.
        //
        // One approach to fixing is to virtually split the spline, which effectively means
        // using p2' = p1 + t(p2 - p1) with t << 1.
        // Another approach is to adjust 's' on the fly to ensure p2 - p0 is not >> ||p1 - p0||
        // We need

        T p0 = *(T*) (p8 + i0 * stride);
        T p1 = *(T*) (p8 + i1 * stride);
        T p2 = *(T*) (p8 + i2 * stride);
        T p3 = *(T*) (p8 + i3 * stride);

        float s = (1.0f - tension) * (1.0f / 6.0f);

        T d1 = s * (p2 - p0);
        T d2 = s * (p3 - p1);

        float l2d1 = sqrlen(d1);
        float l2d2 = sqrlen(d2);

        ratio *= ratio;

        float l1 = sqrlen(p1 - p0) * ratio;
        float l2 = sqrlen(p2 - p1) * ratio;
        float l3 = sqrlen(p3 - p2) * ratio;

        if (l1 > 0)
            l1 = vl_min(l1, l2);
        else
            l1 = l2;

        if (l3 > 0)
            l3 = vl_min(l3, l2);
        else
            l3 = l2;

        if (l2d1 > l1)
            d1 *=  sqrtf(l1 / l2d1);

        if (l2d2 > l3)
            d2 *=  sqrtf(l3 / l2d2);

        T pb1 = p1 + d1;
        T pb2 = p2 - d2;

        return BezierSpline(p1, pb1, pb2, p2);
    }

    template<typename T, typename S> int SplinesFromSamples(int numPoints, const T pi[], S splines[], float tension, size_t stride)
    {
        SL_ASSERT(numPoints >= 0);

        const char* p8 = (const char*) pi;

        switch (numPoints)
        {
        case 0:
            return 0;
        case 1:
            *splines = SplineFromPointsT<T, S>(p8, stride, 0, 0, 0, 0, tension);
            return 1;
        case 2:
            *splines = SplineFromPointsT<T, S>(p8, stride, 0, 0, 1, 1, tension);
            return 1;
        }

        *splines++ = SplineFromPointsT<T, S>(p8, stride, 0, 0, 1, 2, tension);

        for (int i = 0; i < numPoints - 3; i++)
        {
            *splines++ = SplineFromPointsT<T, S>(p8, stride, 0, 1, 2, 3, tension);
            p8 += stride;
        }

        *splines++ = SplineFromPointsT<T, S>(p8, stride, 0, 1, 2, 2, tension);

        return numPoints - 1;
    }

    template<typename T, typename S> int SplinesFromSamplesDynamic(int numPoints, const T pi[], S splines[], float tension, float ratio, size_t stride)
    {
        SL_ASSERT(numPoints >= 0);

        const char* p8 = (const char*) pi;

        switch (numPoints)
        {
        case 0:
            return 0;
        case 1:
            *splines = SplineFromPointsT<T, S>(p8, stride, 0, 0, 0, 0, tension);
            return 1;
        case 2:
            *splines = SplineFromPointsT<T, S>(p8, stride, 0, 0, 1, 1, tension);
            return 1;
        }

        *splines++ = SplineFromPointsDynamicT<T, S>(p8, stride, 0, 0, 1, 2, tension, ratio);

        for (int i = 0; i < numPoints - 3; i++)
        {
            *splines++ = SplineFromPointsDynamicT<T, S>(p8, stride, 0, 1, 2, 3, tension, ratio);
            p8 += stride;
        }

        *splines++ = SplineFromPointsDynamicT<T, S>(p8, stride, 0, 1, 2, 2, tension, ratio);

        return numPoints - 1;
    }

    // Monotone Piecewise Cubic Interpolation
    inline int TriSign(float x, float eps = 1e-6f)
    {
        if (x < -eps)
            return -1;
        if (x > +eps)
            return +1;

        return 0;
    }

    // Given dv = Bezier delta from corresponding adjacent endpoint, and the
    // corresponding end and hull values, return a clamped version of dv that avoids overshoots
    float CubicMonoDelta(float dv, float d0, float d1)
    {
        if (TriSign(d0) != TriSign(d1))
            return 0.0f;

        return Clamp(dv, vl_max(-fabsf(d0), -fabsf(d1)), vl_min(+fabsf(d0), +fabsf(d1)));
    }

    // Reminder that Spline1 holds Bezier coefficients, so
    //   .x = left endpoint
    //   .w = right endpoint
    inline Spline1 MonotonicLeft(float x0, Spline1 s1)
    {
        float dv0 = s1.y - s1.x;

        float d0 = s1.x - x0;
        float d1 = s1.w - s1.x;

        dv0 = CubicMonoDelta(dv0, d0, d1);

        s1.y = s1.x + dv0;
        return s1;
    }

    inline Spline1 MonotonicRight(Spline1 s1, float x2)
    {
        float dv1 = s1.w - s1.z;

        float d1 = s1.w - s1.x;
        float d2 = x2   - s1.w;

        dv1 = CubicMonoDelta(dv1, d1, d2);

        s1.z = s1.w - dv1;
        return s1;
    }

    inline Spline1 Monotonic(float x0, Spline1 s1, float x2)
    {
        float dv0 = s1.y - s1.x;
        float dv1 = s1.w - s1.z;

        float d0 = s1.x - x0;
        float d1 = s1.w - s1.x;
        float d2 = x2   - s1.w;

        dv0 = CubicMonoDelta(dv0, d0, d1);
        dv1 = CubicMonoDelta(dv1, d1, d2);

        s1.y = s1.x + dv0;
        s1.z = s1.w - dv1;
        return s1;
    }

    bool AdvanceAgent(int* indexInOut, float* tInOut, int numSplines)
    // Update index for t if necessary, but don't run off array
    {
        int& index = *indexInOut;
        float& t = *tInOut;

        while (t < 0.0f)
        {
            if (index <= 0)
                return false;

            t += 1.0f;
            index--;
        }

        while (t > 1.0f)
        {
            if (index >= numSplines - 1)
                return false;

            t -= 1.0f;
            index++;
        }

        return true;
    }
}


////////////////////////////////////////////////////////////////////////////////
// 1D
////////////////////////////////////////////////////////////////////////////////


Spline1 SL::BezierSpline(float p0, float p1, float p2, float p3)
{
    return Spline1(p0, p1, p2, p3);
}

Spline1 SL::HermiteSpline(float p0, float p1, float v0, float v1)
    {
    float pb1 = p0 + (1.0f / 3.0f) * v0;
    float pb2 = p1 - (1.0f / 3.0f) * v1;

    return Spline1(p0, pb1, pb2, p1);
}

Spline1 SL::CatmullRomSpline(float p0, float p1, float p2, float p3)
{
    float pb1 = p1 + (1.0f / 6.0f) * (p2 - p0);
    float pb2 = p2 - (1.0f / 6.0f) * (p3 - p1);

    return Spline1(p1, pb1, pb2, p2);
}

Spline1 SL::InterpolatingSpline(float p0, float p1, float p2, float p3)
{
    float pb1 = (1.0f / 3.0f) * p3 - (3.0f / 2.0f) * p2 + 3.0f * p1 - (5.0f / 6.0f) * p0;
    float pb2 = (1.0f / 3.0f) * p0 - (3.0f / 2.0f) * p1 + 3.0f * p2 - (5.0f / 6.0f) * p3;

    return Spline1(p0, pb1, pb2, p3);
}

Spline1 SL::LineSpline(float p0, float p1)
{
    float v = (p1 - p0);

    float pb1 = p0 + (1.0f / 3.0f) * v;
    float pb2 = p1 - (1.0f / 3.0f) * v;

    return Spline1(p0, pb1, pb2, p1);
}

Spline1 SL::CubicSpline(const Vec4f& c)
//Returns Bezier weights for the given cubic coeffs
{
    return Spline1
    (
        c.x,
        c.x + (1.0f / 3.0f) * c.y,
        c.x + (2.0f / 3.0f) * c.y + (1.0f / 3.0f) * c.z,
        c.x +                 c.y +                 c.z + c.w
    );
}

int SL::SplinesFromPoints(int numPoints, const float pi[], Spline1 splines[], float tension, size_t stride)
{
    return ::SplinesFromSamples(numPoints, pi, splines, tension, stride);
}

int SL::SplinesFromPointsDynamic(int numPoints, const float pi[], Spline1 splines[], float tension, float ratio, size_t stride)
{
    return ::SplinesFromSamplesDynamic(numPoints, pi, splines, tension, ratio, stride);
}

int SL::SplinesFromBezier(int numPoints, const float points[], const float hullPoints[], Spline1 splines[], bool split)
{
    int numSplines = split ? numPoints / 2 : numPoints - 1;
    int advance    = split ? 2 : 1;

    for (int i = 0; i < numSplines; i++)
    {
        splines[i] = BezierSpline(points[0], hullPoints[0], hullPoints[1], points[1]);
        points     += advance;
        hullPoints += advance;
    }

    return numSplines;
}

int SL::SplinesFromHermite(int numPoints, const float points[], const float tangents  [], Spline1 splines[], bool split)
{
    int numSplines = split ? numPoints / 2 : numPoints - 1;
    int advance    = split ? 2 : 1;

    for (int i = 0; i < numSplines; i++)
    {
        splines[i] = HermiteSpline(points[0], points[1], tangents[0], tangents[1]);
        points   += advance;
        tangents += advance;
    }

    return numSplines;
}

void SL::MakeMonotonic(int n, Spline1 splines[], bool closed)
{
    if (n >= 2)
    {
        if (closed)
            splines[0] = Monotonic(splines[n - 1].x, splines[0], splines[1].w);
        else
            splines[0] = MonotonicRight(splines[0], splines[1].w);
    }
    else if (n >= 1)
        splines[0] = splines[0];

    for (int i = 1; i < n - 1; i++)
        splines[i] = Monotonic(splines[i - 1].x, splines[i], splines[i + 1].w);

    if (n >= 2)
    {
        if (closed)
            splines[n - 1] = Monotonic(splines[n - 2].x, splines[n - 1], splines[0].w);
        else
            splines[n - 1] = MonotonicLeft(splines[n - 2].x, splines[n - 1]);
    }
}

Vec2f SL::FastBounds(const Spline1& s)
//Returns bounds of the convex hull
{
    Vec2f b01;

    if (s.x <= s.y)
        b01 = Vec2f(s.x, s.y);
    else
        b01 = Vec2f(s.y, s.x);

    Vec2f b23;

    if (s.z <= s.w)
        b23 = Vec2f(s.z, s.w);
    else
        b23 = Vec2f(s.w, s.z);

    return Bounds1
    (
        vl_min(b01.x, b23.x),
        vl_max(b01.y, b23.y)
    );
}

Bounds1 SL::ExactBounds(const Spline1& spline)
//Returns accurate bounds taking extrema into account.
{
    Bounds1 bounds;

    // First take endpoints into account
    if (spline.x <= spline.w)
    {
        bounds.x = spline.x;
        bounds.y = spline.w;
    }
    else
    {
        bounds.x = spline.w;
        bounds.y = spline.x;
    }

    // Now find extrema via standard quadratic equation: c.t' = 0
    Vec4f c = CubicCoeffs(spline);

    float c33 = 3.0f * c.w;
    float cx2 = c.z * c.z - c33 * c.y;

    if (cx2 < 0.0f)
        return bounds;  // no roots!

    float invC33 = 1.0f / c33;
    float ct = -c.z * invC33;
    float cx = sqrtf(cx2) * invC33;

    float t0 = ct + cx;
    float t1 = ct - cx;

    // Must make sure the roots are within the spline interval
    if (t0 > 0.0f && t0 < 1.0f)
    {
        float x = c.x + (c.y + (c.z + c.w * t0) * t0) * t0;

        if      (bounds.x > x)
            bounds.x = x;
        else if (bounds.y < x)
            bounds.y = x;
    }

    if (t1 > 0.0f && t1 < 1.0f)
    {
        float x = c.x + (c.y + (c.z + c.w * t1) * t1) * t1;

        if      (bounds.x > x)
            bounds.x = x;
        else if (bounds.y < x)
            bounds.y = x;
    }

    return bounds;
}

// This is based on one step of De Casteljau's algorithm
void SL::Split(const Spline1& spline, float t, Spline1* spline0, Spline1* spline1)
{
    // assumption: seg = (P0, P1, P2, P3)
    float q0 = lerp(spline.x, spline.y, t);
    float q1 = lerp(spline.y, spline.z, t);
    float q2 = lerp(spline.z, spline.w, t);

    float r0 = lerp(q0, q1, t);
    float r1 = lerp(q1, q2, t);

    float s0 = lerp(r0, r1, t);

    float sx = spline.x;    // support aliasing
    float sw = spline.w;

    *spline0 = Spline1(sx, q0, r0, s0);
    *spline1 = Spline1(s0, r1, q2, sw);
}

// Optimised for t=0.5
void SL::Split(const Spline1& spline, Spline1* spline0, Spline1* spline1)
{
    float q0 = (spline.x + spline.y) * 0.5f;    // x + y / 2
    float q1 = (spline.y + spline.z) * 0.5f;    // y + z / 2
    float q2 = (spline.z + spline.w) * 0.5f;    // z + w / 2

    float r0 = (q0 + q1) * 0.5f;    // x + 2y + z / 4
    float r1 = (q1 + q2) * 0.5f;    // y + 2z + w / 4

    float s0 = (r0 + r1) * 0.5f;    // q0 + 2q1 + q2 / 4 = x+y + 2(y+z) + z+w / 8 = x + 3y + 3z + w

    float sx = spline.x;    // support aliasing
    float sw = spline.w;

    *spline0 = Spline1(sx, q0, r0, s0);
    *spline1 = Spline1(s0, r1, q2, sw);
}

Spline1 SL::Trim(const Spline1& spline_in, float t0, float t1)
{
    SL_ASSERT(t0 <= t1);
    Spline1 spline(spline_in);

    if (t0 > 0.0f)
    {
        float t = t0;

        float q0 = lerp(spline.x, spline.y, t);
        float q1 = lerp(spline.y, spline.z, t);
        float q2 = lerp(spline.z, spline.w, t);

        float r0 = lerp(q0, q1, t);
        float r1 = lerp(q1, q2, t);
        float s0 = lerp(r0, r1, t);

        spline = Spline1(s0, r1, q2, spline.w);
    }

    if (t1 < 1.0f)
    {
        float t = (t1 - t0) / (1.0f - t0);

        float q0 = lerp(spline.x, spline.y, t);
        float q1 = lerp(spline.y, spline.z, t);
        float q2 = lerp(spline.z, spline.w, t);

        float r0 = lerp(q0, q1, t);
        float r1 = lerp(q1, q2, t);
        float s0 = lerp(r0, r1, t);

        spline = Spline1(spline.x, q0, r0, s0);
    }

    return spline;
}

bool SL::Join(const Spline1& s0, const Spline1& s1, Spline1* sOut)
{
    if (s0.w != s1.x) // early out
        return false;

    // assumes t = 0.5

    // backwards solve from left
    float x0 =     s0.x;
    float y0 = 2 * s0.y - x0;
    float z0 = 4 * s0.z - x0 - 2 * y0;
    float w0 = 8 * s0.w - x0 - 3 * (y0 + z0);

    // backwards solve from right
    float w1 =     s1.w;
    float z1 = 2 * s1.z - w1;
    float y1 = 4 * s1.y - w1 - 2 * z1;
    float x1 = 8 * s1.x - w1 - 3 * (y1 + z1);

    float dp = sqr(x0 - x1) + sqr(y0 - y1) + sqr(z0 - z1) + sqr(w0 - w1);

    if (dp < 1e-4f) // do left and right reconstructions agree?
    {
        *sOut = Spline1(x0, y0, z1, w1);   // use most stable terms
        return true;
    }

    return false;
}


////////////////////////////////////////////////////////////////////////////////
// 2D
////////////////////////////////////////////////////////////////////////////////


Spline2 SL::BezierSpline(Vec2f p0, Vec2f p1, Vec2f p2, Vec2f p3)
{
    return Spline2
    {
        Spline1(p0.x, p1.x, p2.x, p3.x),
        Spline1(p0.y, p1.y, p2.y, p3.y),
    };
}

Spline2 SL::HermiteSpline(Vec2f p0, Vec2f p1, Vec2f v0, Vec2f v1)
{
    Vec2f pb1 = p0 + (1.0f / 3.0f) * v0;
    Vec2f pb2 = p1 - (1.0f / 3.0f) * v1;

    return BezierSpline(p0, pb1, pb2, p1);
}

Spline2 SL::CatmullRomSpline(Vec2f p0, Vec2f p1, Vec2f p2, Vec2f p3)
{
    Vec2f pb1 = p1 + (1.0f / 6.0f) * (p2 - p0);
    Vec2f pb2 = p2 - (1.0f / 6.0f) * (p3 - p1);

    return BezierSpline(p1, pb1, pb2, p2);
}

Spline2 SL::InterpolatingSpline(Vec2f p0, Vec2f p1, Vec2f p2, Vec2f p3)
{
    Vec2f pb1 = (1.0f / 3.0f) * p3 - (3.0f / 2.0f) * p2 + 3.0f * p1 - (5.0f / 6.0f) * p0;
    Vec2f pb2 = (1.0f / 3.0f) * p0 - (3.0f / 2.0f) * p1 + 3.0f * p2 - (5.0f / 6.0f) * p3;

    return BezierSpline(p0, pb1, pb2, p3);
}

Spline2 SL::LineSpline(Vec2f p0, Vec2f p1)
{
    Vec2f v = (p1 - p0);

    Vec2f pb1 = p0 + (1.0f / 3.0f) * v;
    Vec2f pb2 = p1 - (1.0f / 3.0f) * v;

    return BezierSpline(p0, pb1, pb2, p1);
}

namespace
{
    const float kCircleOffset = 4.0f / 3.0f * (sqrtf(2.0f) - 1.0f);

    const Spline1 kQuarterB0(1.0f, 1.0f, kCircleOffset, 0.0f);
    const Spline1 kQuarterB1(0.0f, kCircleOffset, 1.0f, 1.0f);
}

Spline2 SL::QuadrantSpline(Vec2f p, float r, int quadrant)
{
    SL_ASSERT(quadrant >= 0 && quadrant < 4);
    Spline2 s;

    switch (quadrant)
    {
    case 0:
        s.xb =  kQuarterB0;
        s.yb =  kQuarterB1;
        break;
    case 1:
        s.xb = -kQuarterB1;
        s.yb =  kQuarterB0;
        break;
    case 2:
        s.xb = -kQuarterB0;
        s.yb = -kQuarterB1;
        break;
    case 3:
        s.xb =  kQuarterB1;
        s.yb = -kQuarterB0;
        break;
    }

    s.xb = s.xb * r + Spline1(p.x, p.x, p.x, p.x);
    s.yb = s.yb * r + Spline1(p.y, p.y, p.y, p.y);

    return s;
}

void SL::CircleSplines(Vec2f p, float r, Spline2 splines[4])
{
    for (int i = 0; i < 4; i++)
        splines[i] = QuadrantSpline(p, r, i);
}

int SL::SplinesFromPoints(int numPoints, const Vec2f pi[], Spline2 splines[], float tension, size_t stride)
{
    return ::SplinesFromSamples(numPoints, pi, splines, tension, stride);
}

int SL::SplinesFromPointsDynamic(int numPoints, const Vec2f pi[], Spline2 splines[], float tension, float ratio, size_t stride)
{
    return ::SplinesFromSamplesDynamic(numPoints, pi, splines, tension, ratio, stride);
}

int SL::SplinesFromBezier(int numPoints, const Vec2f points[], const Vec2f hullPoints[], Spline2 splines[], bool split)
{
    int numSplines = split ? numPoints / 2 : numPoints - 1;
    int advance    = split ? 2 : 1;

    for (int i = 0; i < numSplines; i++)
    {
        splines[i] = BezierSpline(points[0], hullPoints[0], hullPoints[1], points[1]);
        points     += advance;
        hullPoints += advance;
    }

    return numSplines;
}

int SL::SplinesFromHermite(int numPoints, const Vec2f points[], const Vec2f tangents  [], Spline2 splines[], bool split)
{
    int numSplines = split ? numPoints / 2 : numPoints - 1;
    int advance    = split ? 2 : 1;

    for (int i = 0; i < numSplines; i++)
    {
        splines[i] = HermiteSpline(points[0], points[1], tangents[0], tangents[1]);
        points   += advance;
        tangents += advance;
    }

    return numSplines;
}

void SL::MakeMonotonic(int n, Spline2 splines[], bool closed)
{
    if (n >= 2)
    {
        if (closed)
        {
            splines[0].xb = Monotonic(splines[n - 1].xb.x, splines[0].xb, splines[1].xb.w);
            splines[0].yb = Monotonic(splines[n - 1].yb.x, splines[0].yb, splines[1].yb.w);
        }
        else
        {
            splines[0].xb = MonotonicRight(splines[0].xb, splines[1].xb.w);
            splines[0].yb = MonotonicRight(splines[0].yb, splines[1].yb.w);
        }
    }
    else if (n >= 1)
        splines[0] = splines[0];

    for (int i = 1; i < n - 1; i++)
    {
        splines[i].xb = Monotonic(splines[i - 1].xb.x, splines[i].xb, splines[i + 1].xb.w);
        splines[i].yb = Monotonic(splines[i - 1].yb.x, splines[i].yb, splines[i + 1].yb.w);
    }

    if (n >= 2)
    {
        if (closed)
        {
            splines[n - 1].xb = Monotonic(splines[n - 2].xb.x, splines[n - 1].xb, splines[0].xb.w);
            splines[n - 1].yb = Monotonic(splines[n - 2].yb.x, splines[n - 1].yb, splines[0].yb.w);
        }
        else
        {
            splines[n - 1].xb = MonotonicLeft(splines[n - 2].xb.x, splines[n - 1].xb);
            splines[n - 1].yb = MonotonicLeft(splines[n - 2].yb.x, splines[n - 1].yb);
        }
    }
}

namespace
{
    inline Vec2f Evaluate(const Spline2& spline, const Vec4f& w)
    // Evaluate spline with given weights
    {
        return Vec2f
        (
            dot(spline.xb, w),
            dot(spline.yb, w)
        );
    }
}

Vec2f SL::Position(const Spline2& spline, float t)
{
    return Evaluate(spline, BezierWeights(t));
}

Vec2f SL::Velocity(const Spline2& spline, float t)
{
    return Evaluate(spline, BezierWeightsD1(t));
}

Vec2f SL::Acceleration(const Spline2& spline, float t)
{
    return Evaluate(spline, BezierWeightsD2(t));
}

Vec2f SL::Jerk(const Spline2& spline, float)
{
    return Evaluate(spline, kBezierWeightsD3);
}

float SL::Curvature(const Spline2& spline, float t)
{
    Vec2f v = Velocity    (spline, t);
    Vec2f a = Acceleration(spline, t);

    float avCrossLen = fabsf(v.x * a.y - v.y * a.x);
    float vLen = len(v);

    if (vLen == 0.0f)
        return 0.0f;

    return avCrossLen / (vLen * vLen * vLen);
}

Mat2f SL::Frame(const Spline2& spline, float t)
{
    Vec2f v = Velocity    (spline, t);
    Vec2f a = Acceleration(spline, t);

    Mat2f frame;
    frame.x = norm_safe(v);     // x = forward
    frame.y = norm_safe(a - frame.x * dot(frame.x, a));  // y = left or right, pointing in the direction of path turn

    return frame;
}

void SL::Frame(const Spline2& spline, float t, Mat2f* frameInOut)
{
    // Issues:
    //   'b' can flip sign as 'a' moves from one side of the curve to the other
    //   'a' can go to zero on flat sections of the curve. (In which case we should maintain previous values.)
    //   'v' can go to zero at inflection points (ditto).
    // Generally we want to preserve the tangent direction above everything, then the up vector, and last the right vector.
    // E.g., if agent is going into a point and coming back out, we desire that the right vector flips.

    Vec2f v = Velocity    (spline, t);

    Mat2f& frame = *frameInOut;

    Vec2f T;
    if (sqrlen(v) >= 1e-3f)
        T = norm(v);
    else
        T = frame.x;   // fall back to keeping existing tangent

    Vec2f N = cross(T);
    if (dot(N, frame.y) < 0.0f)
        N = -N; // correct b to continue to point same way

    frame.x = T;   // x = forward
    frame.y = N;   // y = left/right
}

Mat2f SL::FrameX(const Spline2& spline, float t)
{
    Vec2f v = Velocity(spline, t);

    Mat2f frame;
    frame.x = norm_safe(v);    // x = forward
    frame.y = cross(frame.x);  // y = left always

    return frame;
}

float SL::LengthEstimate(const Spline2& s, float* error)
{
    // Our convex hull is p0, p1, p2, p3, so p0_p3 is our minimum possible length, and p0_p1 + p1_p2 + p2_p3 our maximum.
    float d03 = sqr(s.xb.x - s.xb.w) + sqr(s.yb.x - s.yb.w);

    float d01 = sqr(s.xb.x - s.xb.y) + sqr(s.yb.x - s.yb.y);
    float d12 = sqr(s.xb.y - s.xb.z) + sqr(s.yb.y - s.yb.z);
    float d23 = sqr(s.xb.z - s.xb.w) + sqr(s.yb.z - s.yb.w);

    float minLength = sqrtf(d03);
    float maxLength = sqrtf(d01) + sqrtf(d12) + sqrtf(d23);

    minLength *= 0.5f;
    maxLength *= 0.5f;

    *error = maxLength - minLength;
    return minLength + maxLength;
}

float SL::Length(const Spline2& s, float maxError)
{
    float error;
    float length = LengthEstimate(s, &error);

    if (error > maxError)
    {
        Spline2 s0;
        Spline2 s1;

        Split(s, &s0, &s1);

        return Length(s0, maxError) + Length(s1, maxError);
    }

    return length;
}

float SL::Length(const Spline2& s, float t0, float t1, float maxError)
{
    SL_ASSERT(t0 >= 0.0f && t0 <= 1.0f);
    SL_ASSERT(t1 >= 0.0f && t1 <= 1.0f);
    SL_ASSERT(t0 <= t1);

    Spline2 s0, s1;

    if (t0 <= 0.0f)
    {
        if (t1 == 1.0f)
            return Length(s, maxError);

        Split(s, t1, &s0, &s1);
        return Length(s0, maxError);
    }
    else if (t1 >= 1.0f)
    {
        Split(s, t0, &s0, &s1);
        return Length(s1, maxError);
    }
    else
    {
        Split(s, t0, &s0, &s1);

        if (t1 == 1.0f)
            return Length(s1, maxError);

        Split(s1, (t1 - t0) / (1.0f - t0), &s0, &s1);
        return Length(s0, maxError);
    }
}

Bounds2 SL::FastBounds(const Spline2& spline)
{
    Vec2f bx = FastBounds(spline.xb);
    Vec2f by = FastBounds(spline.yb);

    return { { bx.x, by.x }, { bx.y, by.y } };
}

Bounds2 SL::ExactBounds(const Spline2& spline)
{
    Vec2f bx = ExactBounds(spline.xb);
    Vec2f by = ExactBounds(spline.yb);

    return { { bx.x, by.x }, { bx.y, by.y } };
}

// Line conversion

int SL::SplinesToLinesLinear(int numSplines, const Spline2 splines[], tEmitLinesFunc2 emitLines, void* context, float step)
{
    Vec2f lineBuffer[kMaxEmitLines][2];
    int lineBufferCount = 0;
    int count = 0;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline2& spline = splines[i];

        lineBuffer[lineBufferCount][0] = Position0(spline);

        for (float t = step; t < 1.0f; t += step)
        {
            lineBuffer[lineBufferCount][1] = Position(spline, t);

            lineBufferCount++;
            if (lineBufferCount == kMaxEmitLines)
            {
                emitLines(lineBufferCount, lineBuffer, context);
                count += lineBufferCount;
                lineBufferCount = 0;
            }
            else
                lineBuffer[lineBufferCount][0] = lineBuffer[lineBufferCount - 1][1];
        }

        lineBuffer[lineBufferCount][1] = Position1(spline);

        lineBufferCount++;
        if (lineBufferCount == kMaxEmitLines)
        {
            emitLines(lineBufferCount, lineBuffer, context);
            count += lineBufferCount;
            lineBufferCount = 0;
        }
    }

    if (lineBufferCount > 0)
    {
        emitLines(lineBufferCount, lineBuffer, context);
        count += lineBufferCount;
    }

    return count;
}

int SL::SplinesToLinesLinear(int numSplines, const Spline2 splines[], tEmitLinesExFunc2 emitLines, void* context, float step)
{
    Vec2f lineBuffer[kMaxEmitLines][2];
    float tvalBuffer[kMaxEmitLines][2];
    int lineBufferCount = 0;
    int count = 0;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline2& spline = splines[i];

        lineBuffer[lineBufferCount][0] = Position0(spline);
        tvalBuffer[lineBufferCount][0] = 0.0f;

        for (float t = step; t < 1.0f; t += step)
        {
            lineBuffer[lineBufferCount][1] = Position(spline, t);
            tvalBuffer[lineBufferCount][1] = t;

            int lastLineBufferCount = lineBufferCount++;

            if (lineBufferCount == kMaxEmitLines)
            {
                emitLines(lineBufferCount, lineBuffer, spline, i, numSplines, tvalBuffer, context);
                count += lineBufferCount;
                lineBufferCount = 0;
            }

            lineBuffer[lineBufferCount][0] = lineBuffer[lastLineBufferCount][1];
            tvalBuffer[lineBufferCount][0] = tvalBuffer[lastLineBufferCount][1];
        }

        lineBuffer[lineBufferCount][1] = Position1(spline);
        tvalBuffer[lineBufferCount][1] = 1.0f;

        lineBufferCount++;
        emitLines(lineBufferCount, lineBuffer, spline, i, numSplines, tvalBuffer, context);
        count += lineBufferCount;
        lineBufferCount = 0;
    }

    return count;
}

namespace
{
    float LineError(const Spline2& s)  // Calculate line approx error as max distance of hull points p1/p2 from p0_p1
    {
        Vec2f p0 = { s.xb.x, s.yb.x };
        Vec2f p1 = { s.xb.y, s.yb.y };
        Vec2f p2 = { s.xb.z, s.yb.z };
        Vec2f p3 = { s.xb.w, s.yb.w };

        Vec2f w = p3 - p0;
        float dww = dot(w, w);

        Vec2f v1 = p1 - p0;
        float dvw1 = dot(v1, w);
        Vec2f v2 = p2 - p0;
        float dvw2 = dot(v2, w);

        float error;

        if (dvw1 <= 0.0f)
            error = sqrlen(v1);
        else if (dvw1 >= dww)
            error = sqrlen(p1 - p3);
        else
            error = sqrlen((dvw1 / dww) * w - v1);

        if (dvw2 <= 0.0f)
            return vl_max(error, sqrlen(v2));
        else if (dvw2 >= dww)
            return vl_max(error, sqrlen(p2 - p3));
        else
            return vl_max(error, sqrlen((dvw2 / dww) * w - v2));
    }
}

int SL::SplinesToLinesAdaptive(int numSplines, const Spline2 splines[], tEmitLinesFunc2 emitLines, void* context, float tolerance)
{
    constexpr int kStackSize = 10; // gives us at least 1024 subdivs -- desirable to cap beyond this anyway.

    Spline2 stack [kStackSize];
    int      stackTop = 0;

    Vec2f lineBuffer[kMaxEmitLines][2];
    int   lineBufferCount = 0;

    int count = 0;
    tolerance *= tolerance;

    for (int i = 0; i < numSplines; i++)
    {
        stackTop = 0;
        stack [stackTop] = splines[i];

        while (stackTop >= 0)
        {
            SL_ASSERT_INDEX(stackTop, kStackSize);
            float error;

            if (stackTop < kStackSize - 1)
                error = LineError(stack[stackTop]);
            else
                error = 0.0f;

            if (error <= tolerance)
            {
                lineBuffer[lineBufferCount][0] = Position0(stack[stackTop]);
                lineBuffer[lineBufferCount][1] = Position1(stack[stackTop]);

                lineBufferCount++;
                if (lineBufferCount == kMaxEmitLines)
                {
                    emitLines(lineBufferCount, lineBuffer, context);
                    lineBufferCount = 0;
                }

                count++;
                stackTop--;
                continue;
            }

            Split(stack[stackTop], stack + stackTop + 1, stack + stackTop + 0);
            stackTop++;
        }
    }

    if (lineBufferCount > 0)
        emitLines(lineBufferCount, lineBuffer, context);

    return count;
}

int SL::SplinesToLinesAdaptive(int numSplines, const Spline2 splines[], tEmitLinesExFunc2 emitLines, void* context, float tolerance)
{
    constexpr int kStackSize = 10; // gives us at least 1024 subdivs -- desirable to cap beyond this anyway.

    Spline2  stack [kStackSize];
    Vec2f    stackT[kStackSize];   // for start/end t values
    int      stackTop = 0;

    Vec2f lineBuffer[kMaxEmitLines][2];
    float tvalBuffer[kMaxEmitLines][2];
    int   lineBufferCount = 0;

    int count = 0;
    tolerance *= tolerance;

    for (int i = 0; i < numSplines; i++)
    {
        stackTop = 0;
        stack [stackTop] = splines[i];
        stackT[stackTop] = { 0.0f, 1.0f };

        while (stackTop >= 0)
        {
            SL_ASSERT_INDEX(stackTop, kStackSize);
            float error;

            if (stackTop < kStackSize - 1)
                error = LineError(stack[stackTop]);
            else
                error = 0.0f;

            if (error <= tolerance)
            {
                lineBuffer[lineBufferCount][0] = Position0(stack[stackTop]);
                lineBuffer[lineBufferCount][1] = Position1(stack[stackTop]);

                tvalBuffer[lineBufferCount][0] = stackT[stackTop].x;
                tvalBuffer[lineBufferCount][1] = stackT[stackTop].y;

                lineBufferCount++;
                if (lineBufferCount == kMaxEmitLines)
                {
                    emitLines(lineBufferCount, lineBuffer, splines[i], i, numSplines, tvalBuffer, context);
                    lineBufferCount = 0;
                }

                count++;
                stackTop--;
                continue;
            }

            Split(stack[stackTop], stack + stackTop + 1, stack + stackTop + 0);
            float halfT = 0.5f * (stackT[stackTop].x + stackT[stackTop].y);
            stackT[stackTop + 1].x = stackT[stackTop + 0].x;
            stackT[stackTop + 1].y = halfT;
            stackT[stackTop + 0].x = halfT;

            stackTop++;
        }

        if (lineBufferCount > 0)
        {
            emitLines(lineBufferCount, lineBuffer, splines[i], i, numSplines, tvalBuffer, context);
            lineBufferCount = 0;
        }
    }

    return count;
}

// Subdivision

void SL::Split(const Spline2& spline, Spline2* spline0, Spline2* spline1)
{
    Split(spline.xb, &spline0->xb, &spline1->xb);
    Split(spline.yb, &spline0->yb, &spline1->yb);
}

void SL::Split(const Spline2& spline, float t, Spline2* spline0, Spline2* spline1)
{
    Split(spline.xb, t, &spline0->xb, &spline1->xb);
    Split(spline.yb, t, &spline0->yb, &spline1->yb);
}

bool SL::Join(const Spline2& s0, const Spline2& s1, Spline2* splineOut)
{
    return
       Join(s0.xb, s1.xb, &splineOut->xb)
    && Join(s0.yb, s1.yb, &splineOut->yb);
}

void SL::Split(vector<Spline2>* splinesIn)
{
    vector<Spline2> splines;

    for (const Spline2& s : *splinesIn)
    {
        Spline2 s0, s1;

        Split(s, &s0, &s1);
        splines.push_back(s0);
        splines.push_back(s1);
    }

    splinesIn->swap(splines);
}

void SL::Split(vector<Spline2>* splinesIn, int n)
{
    vector<Spline2> splines;

    for (const Spline2& s : *splinesIn)
    {
        Spline2 ss(s);
        Spline2 s0, s1;

        for (int i = n; i > 1; i--)
        {
            Split(ss, 1.0f / i, &s0, &ss);
            splines.push_back(s0);
        }
        splines.push_back(ss);
    }

    splinesIn->swap(splines);
}

void SL::Join(vector<Spline2>* splinesIn)
{
    vector<Spline2> splines;
    const Spline2* prevS = 0;

    for (const Spline2& s : *splinesIn)
    {
        if (!prevS)
        {
            prevS = &s;
            continue;
        }

        Spline2 sj;
        if (Join(*prevS, s, &sj))
            splines.push_back(sj);
        else
        {
            splines.push_back(*prevS);
            splines.push_back(s);
        }

        prevS = 0;
    }

    if (prevS)
        splines.push_back(*prevS);

    splinesIn->swap(splines);
}

Spline2 SL::Trim(const Spline2& spline, float t0, float t1)
{
    return { Trim(spline.xb, t0, t1), Trim(spline.yb, t0, t1) };
}

namespace
{
    void SubdivideForLength(const Spline2& s, vector<Spline2>* splines, float tolerance)
    {
        float error;
        float length = LengthEstimate(s, &error);

        if (error <= tolerance * length)
            splines->push_back(s);
        else
        {
            Spline2 s1, s2;
            Split(s, &s1, &s2);

            SubdivideForLength(s1, splines, tolerance);
            SubdivideForLength(s2, splines, tolerance);
        }
    }
}

void SL::SubdivideForLength(vector<Spline2>* splinesIn, float tolerance)
{
    vector<Spline2> splines;

    for (const Spline2& s : *splinesIn)
        ::SubdivideForLength(s, &splines, tolerance);

    splinesIn->swap(splines);
}

namespace
{
    inline float ArcError(const Spline2& s, float* tSplit)
    {
        Vec2f ex = ArcError2(s.xb);
        Vec2f ey = ArcError2(s.yb);

        float e0 = ex.x + ey.x;
        float e1 = ex.y + ey.y;
        float es2 = e0 + e1;

        float f = (es2 < 1e-6f) ? 0.5f : sqrtf(e0 / es2);

        *tSplit = (1.0f / 3.0f) * (1.0f + f);

        return sqrtf(es2);
    }

    void SubdivideForT(const Spline2& s, vector<Spline2>* splines, float tolerance)
    {
        float splitT;
        float err = ArcError(s, &splitT);

        if (err <= tolerance)
            splines->push_back(s);
        else
        {
            Spline2 s1, s2;
            Split(s, splitT, &s1, &s2);

            SubdivideForT(s1, splines, tolerance);
            SubdivideForT(s2, splines, tolerance);
        }
    }
}

void SL::SubdivideForT(vector<Spline2>* splinesIn, float tolerance)
{
    vector<Spline2> splines;

    for (const Spline2& s : *splinesIn)
        ::SubdivideForT(s, &splines, tolerance);

    splinesIn->swap(splines);
}

void SL::SubdivideForLengthRatio(vector<Vec2f>& positions, float maxRatio)
{
    for (int i = 1; i < size_i(positions); )
    {
        float l1 = len(positions[i] - positions[i - 1]);

        if (i > 1)
        {
            float l0 = len(positions[i - 1] - positions[i - 2]);

            if (l1 > l0 * maxRatio)
            {
                positions.insert(positions.begin() + i, 0.5f * (positions[i] + positions[i - 1]));
                continue;
            }
        }

        if (i < size_i(positions) - 1)
        {
            float l2 = len(positions[i + 1] - positions[i]);

            if (l1 > l2 * maxRatio)
            {
                positions.insert(positions.begin() + i, 0.5f * (positions[i] + positions[i - 1]));
                continue;
            }
        }

        i++;
    }
}

namespace
{
    inline float ClosestPoint(Vec2f p, Vec2f p0, Vec2f p1)
    {
        Vec2f w = p1 - p0;
        Vec2f v =  p - p0;

        float dvw = dot(v, w);

        if (dvw <= 0.0f)
            return 0.0f;

        float dww = dot(w, w);

        if (dvw >= dww)
            return 1.0f;

        return dvw / dww;
    }

    void FindClosestPointNewtonRaphson(const Spline2& spline, Vec2f p, float sIn, int maxIterations, float* tOut, float* dOut)
    {
        SL_ASSERT(sIn >= 0.0f && sIn <= 1.0f);
        const float maxS = 1.0f - 1e-6f;

        float sk = sIn;
        float dk = len(Position(spline, sk) - p);

        constexpr float width = 1e-3f;

        float maxJump  = 0.5f;   // avoid jumping too far, leads to oscillation

        for (int i = 0; i < maxIterations; i++)
        {
            float ss = Clamp(sk, width, 1.0f - width); // so can interpolate points for Newtons method

            float d1 = len(Position(spline, ss - width) - p);
            float d2 = len(Position(spline, ss        ) - p);
            float d3 = len(Position(spline, ss + width) - p);

            float g1 = (d2 - d1) / width;
            float g2 = (d3 - d2) / width;

            float grad = (d3 - d1) / (2.0f * width);
            float curv = (g2 - g1) / width;

            float sn;   // next candidate

            if (curv > 0.0f)    // if d' is heading towards a minima, apply NR for d'
                sn = ss - grad / curv;
            else if (grad != 0.0f)
                sn = ss - d2 / grad; // otherwise, apply for D.
            else
                sn = sk;

            sn = Clamp(sn, sk - maxJump, sk + maxJump);   // avoid large steps, often unstable.

            // only update our estimate if the new value is in range and closer.
            if (sn >= 0.0f && sn < maxS)
            {
                float dn = len(Position(spline, sn) - p);

                if (dn < dk)    // only update sk if d actually gets smaller
                {
                    sk = sn;
                    dk = dn;
                }
            }

            maxJump *= 0.5f;    // reduce on a schedule -- helps binary search back towards a jump that is valid.
        }

        (*tOut) = sk;
        (*dOut) = dk;
    }
}

float SL::FindClosestPoint(Vec2f p, const Spline2& spline)
{
    // Approximate s from straight line between the start and end.
    float s = ClosestPoint(p, Position0(spline), Position1(spline));

    // Use this as starting point for Newton-Raphson solve.
    float d;
    FindClosestPointNewtonRaphson(spline, p, s, 8, &s, &d);

    return s;
}

float SL::FindClosestPoint(Vec2f p, int numSplines, const Spline2 splines[], int* index)
{
    vector<SubSpline2> nearbyInfo;

    FindNearbySplines(p, numSplines, splines, &nearbyInfo);
    return FindClosestPoint(p, numSplines, splines, nearbyInfo, index);
}


namespace
{
    void FindMinMaxDistance2s(Vec2f p, const Bounds2& bbox, float* minD2, float* maxD2)
    {
        const Vec2f& p0 = bbox.mMin;
        const Vec2f& p1 = bbox.mMax;

        // Find the nearest point to p inside the bbox
        // This can be a bbox vertex, a point on an edge or face, or p itself if it's inside the box
        float minX = Clamp(p.x, p0.x, p1.x);
        float minY = Clamp(p.y, p0.y, p1.y);

        // Find the farthest point from p inside the bbox
        // This is always a bbox vertex.
        Vec2f d0(abs(p - p0));
        Vec2f d1(abs(p - p1));

        float maxX = d0.x > d1.x ? p0.x : p1.x; // select the coordinate we're farthest from
        float maxY = d0.y > d1.y ? p0.y : p1.y;

        // return len2
        *minD2 = sqr(p.x - minX) + sqr(p.y - minY);
        *maxD2 = sqr(p.x - maxX) + sqr(p.y - maxY);
    }

    void FindMinMaxDistance2s(Vec2f p, const Spline2& spline, float* minD2, float* maxD2)
    {
        Bounds2 bbox = FastBounds(spline);

        FindMinMaxDistance2s(p, bbox, minD2, maxD2);
    }

    void Split(const SubSpline2& s, SubSpline2* s0, SubSpline2* s1)
    {
        ::Split(s.mSpline.xb, &s0->mSpline.xb, &s1->mSpline.xb);
        ::Split(s.mSpline.yb, &s0->mSpline.yb, &s1->mSpline.yb);

        s0->mParent = s.mParent;
        s1->mParent = s.mParent;
    }
}

int SL::FindNearbySplines(Vec2f p, int numSplines, const Spline2 splines[], vector<SubSpline2>* results, float* smallestFarOut, int numIter)
{
    vector<SubSpline2>& nearSplines = *results;

    nearSplines.clear();

    float smallestFar  = FLT_MAX;
    float smallestNear = FLT_MAX;

    // Find initial list

    int maxSize = 0;

    for (int i = 0; i < numSplines; i++)
    {
        float near;
        float far;
        FindMinMaxDistance2s(p, splines[i], &near, &far);

        if (near < smallestFar)
        {
            // we at least overlap the current set.
            if (near < smallestNear)
                smallestNear = near;

            if (far < smallestFar)
            {
                // we bring in the 'best' far distance
                smallestFar = far;

                // compact list to reject any segments that now cannot be closest.
                int dj = 0;

                for (int j = 0, nj = size_i(nearSplines); j < nj; j++)
                    if (nearSplines[j].mD2 <= smallestFar)
                    {
                        if (dj < j)
                            nearSplines[dj] = nearSplines[j];

                        dj++;
                    }

                nearSplines.resize(dj);
            }

            SubSpline2 ss = { splines[i], i, near };
            nearSplines.push_back(ss);

            if (maxSize < size_i(nearSplines))
                maxSize = size_i(nearSplines);
        }
    }

    // Subdivide + refine

    int numNearSplines = size_i(nearSplines);

    for (int i = 0; i < numIter; i++)
    {
        int numNearSplines2 = numNearSplines * 2;

        nearSplines.resize(numNearSplines2);

        for (int i = numNearSplines - 1; i >= 0; i--)
            ::Split(nearSplines[i], &nearSplines[2 * i], &nearSplines[2 * i + 1]);

        smallestNear = FLT_MAX; // this may actually increase on subdivision.

        for (int i = 0; i < numNearSplines2; i++)
        {
            float near;
            float far;
            FindMinMaxDistance2s(p, nearSplines[i].mSpline, &near, &far);

            if (far < smallestFar)
                smallestFar = far;
            if (near < smallestNear)
                smallestNear = near;

            nearSplines[i].mD2 = near;
        }

        int di = 0;
        for (int i = 0; i < numNearSplines2; i++)
            if (nearSplines[i].mD2 < smallestFar)
            {
                if (di < i)
                    nearSplines[di] = nearSplines[i];

                di++;
            }

        nearSplines.resize(di);
        numNearSplines = di;
    }

    SL_ASSERT(numNearSplines > 0 || numSplines == 0);

    if (smallestFarOut)
        *smallestFarOut = smallestFar;

    return numNearSplines;
}

float SL::FindClosestPoint(Vec2f p, int numSplines, const Spline2 splines[], const vector<SubSpline2>& nearbySplines, int* index)
{
    int prevParent = -1;
    float minD = FLT_MAX;
    float minT = 0.0f;

    *index = -1;

    for (const SubSpline2& subSpline : nearbySplines)
    {
        if (subSpline.mParent != prevParent)
        {
            SL_ASSERT(subSpline.mParent >= 0 && subSpline.mParent < numSplines);
            (void) numSplines;

            const Spline2& spline = splines[subSpline.mParent];

            float t = ClosestPoint(p, Position0(spline), Position1(spline));
            float d;

            FindClosestPointNewtonRaphson(spline, p, t, 8, &t, &d);

            if (minD > d)
            {
                minD = d;
                *index = subSpline.mParent;
                minT = t;
            }
        }
    }

    return minT;
}

namespace
{
    struct SubSplineT2
    {
        Spline2 mSpline;
        float    mT0;
        float    mT1;
    };

    inline void Split(const SubSplineT2& s, SubSplineT2* s0, SubSplineT2* s1)
    {
        ::Split(s.mSpline.xb, &s0->mSpline.xb, &s1->mSpline.xb);
        ::Split(s.mSpline.yb, &s0->mSpline.yb, &s1->mSpline.yb);

        float midT = 0.5f * (s.mT0 + s.mT1);

        s0->mT0     = s.mT0;
        s0->mT1     = midT;

        s1->mT0     = midT;
        s1->mT1     = s.mT1;
    }

    int FindSubSplineIntersections
    (
        const SubSplineT2& spline0,
        const SubSplineT2& spline1,
        int   dest,
        int   maxDest,
        float results[][2],
        float tolerance
    )
    {
        SL_ASSERT(dest < maxDest);

        Bounds2 bbox0 = ExactBounds(spline0.mSpline);
        Bounds2 bbox1 = ExactBounds(spline1.mSpline);

        if (!Intersects(bbox0, bbox1))
            return dest;

        if (Larger(bbox0, tolerance))
        {
            SubSplineT2 spline00, spline01;
            Split(spline0, &spline00, &spline01);

            if (Larger(bbox1, tolerance))
            {
                SubSplineT2 spline10, spline11;
                Split(spline1, &spline10, &spline11);

                dest     = FindSubSplineIntersections(spline00, spline10, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline01, spline10, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline00, spline11, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline01, spline11, dest, maxDest, results, tolerance);
            }
            else
            {
                dest     = FindSubSplineIntersections(spline00, spline1, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline01, spline1, dest, maxDest, results, tolerance);
            }

            return dest;
        }

        if (Larger(bbox1, tolerance))
        {
            SubSplineT2 spline10, spline11;
            Split(spline1, &spline10, &spline11);

            dest     = FindSubSplineIntersections(spline0, spline10, dest, maxDest, results, tolerance);
            if (dest < maxDest)
                dest = FindSubSplineIntersections(spline0, spline11, dest, maxDest, results, tolerance);

            return dest;
        }

        float t0 = 0.5f * (spline0.mT0 + spline0.mT1);
        float t1 = 0.5f * (spline1.mT0 + spline1.mT1);

        // debounce
        for (int i = 0; i < dest; i++)
            if (fabsf(results[i][0] - t0) < 1e-2f
             && fabsf(results[i][1] - t1) < 1e-2f)
            {
                return dest;
            }

        results[dest][0] = t0;
        results[dest][1] = t1;
        return dest + 1;
    }
}

int SL::FindSplineIntersections(const Spline2& spline0, const Spline2& spline1, int maxResults, float results[][2], float tolerance)
{
    if (maxResults <= 0)
        return 0;

    SubSplineT2 subSpline0 = { spline0, 0.0f, 1.0f };
    SubSplineT2 subSpline1 = { spline1, 0.0f, 1.0f };

    return FindSubSplineIntersections(subSpline0, subSpline1, 0, maxResults, results, tolerance);
}

namespace
{
    int FindSplineIntersections
    (
        int   is0, int numSplines0, const Spline2 splines0[],
        int   is1, int numSplines1, const Spline2 splines1[],
        int   maxResults,
        int   resultsI[][2],
        float resultsT[][2],
        float tolerance
    )
    {
        if (maxResults <= 0 || numSplines0 == 0 || numSplines1 == 0)
            return 0;

        int count = 0;

        // once the lists are small enough, brute-force the remainder, as recalculating the bounds is not free
        if (numSplines0 >= 4 || numSplines1 >= 4)
        {
            // Terminate if the lists don't overlap
            Bounds2 b0 = FastBounds(splines0[0]);
            Bounds2 b1 = FastBounds(splines1[0]);

            for (int i = 1; i < numSplines0; i++)
                Add(b0, FastBounds(splines0[i]));
            for (int i = 1; i < numSplines1; i++)
                Add(b1, FastBounds(splines1[i]));

            if (!Intersects(b0, b1))
                return 0;

            // Divide each spline list into two, and recurse
            int n00 = numSplines0 / 2;
            int n01 = numSplines0 - n00;
            int n10 = numSplines1 / 2;
            int n11 = numSplines1 - n10;

            count += FindSplineIntersections(is0 +   0, n00, splines0 +   0, is1 +   0, n10, splines1 +   0, maxResults - count, resultsI + count, resultsT + count, tolerance);
            count += FindSplineIntersections(is0 +   0, n00, splines0 +   0, is1 + n10, n11, splines1 + n10, maxResults - count, resultsI + count, resultsT + count, tolerance);
            count += FindSplineIntersections(is0 + n00, n01, splines0 + n00, is1 +   0, n10, splines1 +   0, maxResults - count, resultsI + count, resultsT + count, tolerance);
            count += FindSplineIntersections(is0 + n00, n01, splines0 + n00, is1 + n10, n11, splines1 + n10, maxResults - count, resultsI + count, resultsT + count, tolerance);

            return count;
        }

        SubSplineT2 st0, st1;

        st0.mT0 = 0.0f;
        st0.mT1 = 1.0f;
        st1.mT0 = 0.0f;
        st1.mT1 = 1.0f;

        for (int i0 = 0; i0 < numSplines0; i0++)
            for (int i1 = 0; i1 < numSplines1; i1++)
            {
                st0.mSpline = splines0[i0];
                st1.mSpline = splines1[i1];

                int numIntersections = FindSubSplineIntersections(st0, st1, 0, maxResults - count, resultsT + count, tolerance);

                for (int k = 0; k < numIntersections; k++)
                {
                    resultsI[k + count][0] = is0 + i0;
                    resultsI[k + count][1] = is1 + i1;
                }

                count += numIntersections;
                SL_ASSERT(count <= maxResults);

                if (count == maxResults)
                    return count;
            }

        return count;
    }
}

int SL::FindSplineIntersections
(
    int   numSplines0, const Spline2 splines0[],
    int   numSplines1, const Spline2 splines1[],
    int   maxResults,
    int   resultsI[][2],
    float resultsT[][2],
    float tolerance
)
{
    return ::FindSplineIntersections(0, numSplines0, splines0, 0, numSplines1, splines1, maxResults, resultsI, resultsT, tolerance);
}

int SL::FindSplineIntersections
(
    int   numSplines, const Spline2 splines[],
    int   maxResults,
    int   resultsI[][2],
    float resultsT[][2],
    float tolerance
)
// Find self-intersections
{
    if (maxResults <= 0 || numSplines == 0)
        return 0;

    int count = 0;

    if (numSplines >= 8)
    {
        const int n0 = numSplines / 2;
        const int n1 = numSplines - n0;

        // Find intersections between the two halves
        count += ::FindSplineIntersections(0, n0, splines, 0, n1, splines + n0, maxResults - count, resultsI + count, resultsT + count, tolerance);

        for (int i = 0; i < count; i++)
        {
            // ignore spurious intersections between endpoint of first and start of second
            if (resultsI[i][1] == 0 && resultsI[i][0] == n0 - 1)
                if (resultsT[i][0] > 0.95f && resultsT[i][1] < 0.05f)
                {
                    resultsT[i][0] = resultsT[count - 1][0];
                    resultsT[i][1] = resultsT[count - 1][1];
                    resultsI[i][0] = resultsI[count - 1][0];
                    resultsI[i][1] = resultsI[count - 1][1];

                    i--;
                    count--;
                    continue;
                }

            resultsI[i][1] += n0;           // adjust second half indices
        }

        // Find self-intersections in the first half
        count += FindSplineIntersections(n0, splines +  0, maxResults - count, resultsI + count, resultsT + count, tolerance);

        // Find self-intersections in the second half
        int prevCount = count;
        count += FindSplineIntersections(n1, splines + n0, maxResults - count, resultsI + count, resultsT + count, tolerance);

        for (int i = prevCount; i < count; i++)
        {
            resultsI[i][0] += n0;   // adjust second half indices
            resultsI[i][1] += n0;
        }

        return count;
    }

    SubSplineT2 st0, st1;

    st0.mT0 = 0.0f;
    st0.mT1 = 1.0f;
    st1.mT0 = 0.0f;
    st1.mT1 = 1.0f;

    for (int i0 = 0; i0 < numSplines; i0++)
        for (int i1 = i0 + 1; i1 < numSplines; i1++)
        {
            st0.mSpline = splines[i0];
            st1.mSpline = splines[i1];

            int numIntersections = FindSubSplineIntersections(st0, st1, 0, maxResults, resultsT, tolerance);

            if (i1 == i0 + 1)   // ignore spurious intersections between endpoint of i0 and start of i1
            {
                for (int k = 0; k < numIntersections; k++)
                {
                    if (resultsT[k][0] > 0.95f && resultsT[k][1] < 0.05f)
                    {
                        resultsT[k][0] = resultsT[numIntersections - 1][0];
                        resultsT[k][1] = resultsT[numIntersections - 1][1];

                        k--;
                        numIntersections--;
                    }
                }
            }

            for (int k = 0; k < numIntersections; k++)
            {
                resultsI[k][0] = i0;
                resultsI[k][1] = i1;
            }

            count += numIntersections;
            maxResults -= numIntersections;

            if (maxResults <= 0)
                return count;

            resultsT += numIntersections;
            resultsI += numIntersections;
        }

    return count;
}

float SL::AdvanceAgent(const Spline2& spline, float t, float ds)
{
    Vec2f v  = Velocity(spline, t);
    float v2 = sqrlen(v);
    float dt = ds;

    if (v2 > 0.01f)
        dt *= InvSqrtFast(v2);
    else
        dt *= 10.0f;

    return t + dt;
}

bool SL::AdvanceAgent(int* index, float* t, int numSplines, const Spline2 splines[], float ds)
{
    *t = AdvanceAgent(splines[*index], *t, ds);

    return ::AdvanceAgent(index, t, numSplines);
}

Spline2 SL::Reverse(const Spline2& spline)
{
    return
    {
        Vec4f(spline.xb.w, spline.xb.z, spline.xb.y, spline.xb.x),
        Vec4f(spline.yb.w, spline.yb.z, spline.yb.y, spline.yb.x)
    };
}

void SL::Reverse(vector<Spline2>* splines)
{
    int n = int(splines->size());
    int h = n / 2;

    for (int i = 0; i < h; i++)
    {
        Spline2& s0 = (*splines)[i];
        Spline2& s1 = (*splines)[n - i - 1];

        Spline2 sr0 = Reverse(s1);
        Spline2 sr1 = Reverse(s0);

        s0 = sr0;
        s1 = sr1;
    }

    if (2 * h < n)
        (*splines)[h] = Reverse((*splines)[h]);
}

Spline2 SL::Offset(const Spline2& spline, float offset)
{
    float sx0 = spline.xb.y - spline.xb.x;
    float sy0 = spline.yb.y - spline.yb.x;
    float sd0 = InvSqrtFast(sx0 * sx0 + sy0 * sy0) * offset;

    float sx1 = spline.xb.z - spline.xb.x;
    float sy1 = spline.yb.z - spline.yb.x;
    float sd1 = InvSqrtFast(sx1 * sx1 + sy1 * sy1) * offset;

    float sx2 = spline.xb.w - spline.xb.y;
    float sy2 = spline.yb.w - spline.yb.y;
    float sd2 = InvSqrtFast(sx2 * sx2 + sy2 * sy2) * offset;

    float sx3 = spline.xb.w - spline.xb.z;
    float sy3 = spline.yb.w - spline.yb.z;
    float sd3 = InvSqrtFast(sx3 * sx3 + sy3 * sy3) * offset;

    return
    {
        spline.xb + Vec4f(sy0 * sd0, sy1 * sd1, sy2 * sd2, sy3 * sd3),
        spline.yb - Vec4f(sx0 * sd0, sx1 * sd1, sx2 * sd2, sx3 * sd3)
    };
}

void SL::Offset(vector<Spline2>* splines, float offset)
{
    for (Spline2& s : *splines)
        s = Offset(s, offset);
}


////////////////////////////////////////////////////////////////////////////////
// 3D
////////////////////////////////////////////////////////////////////////////////


Spline3 SL::BezierSpline(Vec3f p0, Vec3f p1, Vec3f p2, Vec3f p3)
{
    return Spline3
    {
        Vec4f(p0.x, p1.x, p2.x, p3.x),
        Vec4f(p0.y, p1.y, p2.y, p3.y),
        Vec4f(p0.z, p1.z, p2.z, p3.z)
    };
}

Spline3 SL::HermiteSpline(Vec3f p0, Vec3f p1, Vec3f v0, Vec3f v1)
{
    Vec3f pb1 = p0 + (1.0f / 3.0f) * v0;
    Vec3f pb2 = p1 - (1.0f / 3.0f) * v1;

    return BezierSpline(p0, pb1, pb2, p1);
}

Spline3 SL::CatmullRomSpline(Vec3f p0, Vec3f p1, Vec3f p2, Vec3f p3)
{
    Vec3f pb1 = p1 + (1.0f / 6.0f) * (p2 - p0);
    Vec3f pb2 = p2 - (1.0f / 6.0f) * (p3 - p1);

    return BezierSpline(p1, pb1, pb2, p2);
}

Spline3 SL::InterpolatingSpline(Vec3f p0, Vec3f p1, Vec3f p2, Vec3f p3)
{
    Vec3f pb1 = (1.0f / 3.0f) * p3 - (3.0f / 2.0f) * p2 + 3.0f * p1 - (5.0f / 6.0f) * p0;
    Vec3f pb2 = (1.0f / 3.0f) * p0 - (3.0f / 2.0f) * p1 + 3.0f * p2 - (5.0f / 6.0f) * p3;

    return BezierSpline(p0, pb1, pb2, p3);
}

Spline3 SL::LineSpline(Vec3f p0, Vec3f p1)
{
    Vec3f v = (p1 - p0);

    Vec3f pb1 = p0 + (1.0f / 3.0f) * v;
    Vec3f pb2 = p1 - (1.0f / 3.0f) * v;

    return BezierSpline(p0, pb1, pb2, p1);
}

Spline3 SL::QuadrantSpline(Vec3f p, float r, int quadrant)
{
    Spline3 result;
    (Spline2&) result = QuadrantSpline({ p.x, p.y }, r, quadrant);
    result.zb = Spline1(p.z, p.z, p.z, p.z);
    return result;
}

void SL::CircleSplines(Vec3f p, float r, Spline3 splines[4])
{
    for (int i = 0; i < 4; i++)
        splines[i] = QuadrantSpline(p, r, i);
}

int SL::SplinesFromPoints(int numPoints, const Vec3f pi[], Spline3 splines[], float tension, size_t stride)
{
    return ::SplinesFromSamples(numPoints, pi, splines, tension, stride);
}

int SL::SplinesFromPointsDynamic(int numPoints, const Vec3f pi[], Spline3 splines[], float tension, float ratio, size_t stride)
{
    return ::SplinesFromSamplesDynamic(numPoints, pi, splines, tension, ratio, stride);
}

int SL::SplinesFromBezier(int numPoints, const Vec3f points[], const Vec3f hullPoints[], Spline3 splines[], bool split)
{
    int numSplines = split ? numPoints / 2 : numPoints - 1;
    int advance    = split ? 2 : 1;

    for (int i = 0; i < numSplines; i++)
    {
        splines[i] = BezierSpline(points[0], hullPoints[0], hullPoints[1], points[1]);
        points     += advance;
        hullPoints += advance;
    }

    return numSplines;
}

int SL::SplinesFromHermite(int numPoints, const Vec3f points[], const Vec3f tangents  [], Spline3 splines[], bool split)
{
    int numSplines = split ? numPoints / 2 : numPoints - 1;
    int advance    = split ? 2 : 1;

    for (int i = 0; i < numSplines; i++)
    {
        splines[i] = HermiteSpline(points[0], points[1], tangents[0], tangents[1]);
        points   += advance;
        tangents += advance;
    }

    return numSplines;
}

void SL::MakeMonotonic(int n, Spline3 splines[], bool closed)
{
    if (n >= 2)
    {
        if (closed)
        {
            splines[0].xb = Monotonic(splines[n - 1].xb.x, splines[0].xb, splines[1].xb.w);
            splines[0].yb = Monotonic(splines[n - 1].yb.x, splines[0].yb, splines[1].yb.w);
            splines[0].zb = Monotonic(splines[n - 1].zb.x, splines[0].yb, splines[1].zb.w);
        }
        else
        {
            splines[0].xb = MonotonicRight(splines[0].xb, splines[1].xb.w);
            splines[0].yb = MonotonicRight(splines[0].yb, splines[1].yb.w);
            splines[0].zb = MonotonicRight(splines[0].zb, splines[1].zb.w);
        }
    }
    else if (n >= 1)
        splines[0] = splines[0];

    for (int i = 1; i < n - 1; i++)
    {
        splines[i].xb = Monotonic(splines[i - 1].xb.x, splines[i].xb, splines[i + 1].xb.w);
        splines[i].yb = Monotonic(splines[i - 1].yb.x, splines[i].yb, splines[i + 1].yb.w);
        splines[i].zb = Monotonic(splines[i - 1].zb.x, splines[i].zb, splines[i + 1].zb.w);
    }

    if (n >= 2)
    {
        if (closed)
        {
            splines[n - 1].xb = Monotonic(splines[n - 2].xb.x, splines[n - 1].xb, splines[0].xb.w);
            splines[n - 1].yb = Monotonic(splines[n - 2].yb.x, splines[n - 1].yb, splines[0].yb.w);
            splines[n - 1].zb = Monotonic(splines[n - 2].zb.x, splines[n - 1].yb, splines[0].zb.w);
        }
        else
        {
            splines[n - 1].xb = MonotonicLeft(splines[n - 2].xb.x, splines[n - 1].xb);
            splines[n - 1].yb = MonotonicLeft(splines[n - 2].yb.x, splines[n - 1].yb);
            splines[n - 1].zb = MonotonicLeft(splines[n - 2].zb.x, splines[n - 1].zb);
        }
    }
}

namespace
{
    inline Vec3f Evaluate(const Spline3& spline, const Vec4f& w)
    // Evaluate spline with given weights
    {
        return Vec3f
        (
            dot(spline.xb, w),
            dot(spline.yb, w),
            dot(spline.zb, w)
        );
    }
}

Vec3f SL::Position(const Spline3& spline, float t)
{
    return Evaluate(spline, BezierWeights(t));
}

Vec3f SL::Velocity(const Spline3& spline, float t)
{
    return Evaluate(spline, BezierWeightsD1(t));
}

Vec3f SL::Acceleration(const Spline3& spline, float t)
{
    return Evaluate(spline, BezierWeightsD2(t));
}

Vec3f SL::Jerk(const Spline3& spline, float)
{
    return Evaluate(spline, kBezierWeightsD3);
}

float SL::Curvature(const Spline3& spline, float t)
{
    Vec3f v = Velocity    (spline, t);
    Vec3f a = Acceleration(spline, t);

    float avCrossLen = len(cross(v, a));
    float vLen = len(v);

    if (vLen == 0.0f)
        return 0.0f;

    return avCrossLen / (vLen * vLen * vLen);
}

float SL::Torsion(const Spline3& spline, float t)
{
    Vec3f v = Velocity    (spline, t);
    Vec3f a = Acceleration(spline, t);
    Vec3f j = Jerk        (spline, t);

    float avCrossLen2 = sqrlen(cross(v, a));

    if (avCrossLen2 == 0.0f)
        return 0.0f;

    return dot(cross(v, a), j) / avCrossLen2;
}

Mat3f SL::Frame(const Spline3& spline, float t)
{
    Vec3f v = Velocity    (spline, t);  // r'
    Vec3f a = Acceleration(spline, t);  // r''

    Vec3f T = norm_safe(v);             // tangent
    Vec3f B = norm_safe(cross(v, a));   // binormal
    Vec3f N = cross(B, T);              // reconstructed 'normal' orthogonal to tangent

    Mat3f frame;
    frame.x = T;   // x = forward
    frame.y = N;   // y = left/right, direction of curvature, defines plane of curvature
    frame.z = B;   // z = up/down, perpendicular to plane of curvature, direction of torsion

    return frame;
}

void SL::Frame(const Spline3& spline, float t, Mat3f* frameInOut)
{
    // The goal here is to move the incoming frame the least amount possible to align
    // with what T/N/B are telling us. While x always has to match the tangent, we
    // allow ambiguity between y and z, so e.g. an upwards curving path doesn't
    // immediately flip z-up to to z-right.

    // Issues:
    //   'B' can flip direction as 'a' moves from one side of the curve to the other
    //   'a' can go to zero on flat sections of the curve. (In which case we should maintain previous values.)
    //   'v' can go to zero at inflection points (ditto).
    // Generally we want to preserve the tangent direction above everything, then the up vector, and last the right vector.
    // E.g., if agent is going into a point and coming back out, we desire that the right vector flips rather than up.

    Vec3f v = Velocity    (spline, t);
    Vec3f a = Acceleration(spline, t);

    Mat3f& frame = *frameInOut;

    Vec3f T;
    if (sqrlen(v) >= 1e-3f)
        T = norm(v);
    else
        T = frame.x;   // fall back to keeping existing tangent

    Vec3f ga = cross(T, a);
    Vec3f B, N;

    if (sqrlen(ga) >= 1e-3f)
    {
        B = norm(ga);

        float dp = dot(B, frame.z);

        if (fabsf(dp) >= 0.7f)  // ~cos45, we're playing a bit fast and loose as we expect the frame to evolve slowly...
        {
            if (dp < 0.0f)
                B = -B; // correct B to continue to point same way

            N = cross(B, T);
        }
        else
        {
            N = B;
            B = cross(T, N);

            dp = dot(B, frame.z);  // still prioritize B pointing the same way...

            if (dp < 0.0f)
            {
                N = -N;
                B = -B;
            }
        }
    }
    else
    {
        B = frame.z;   // fall back to previous up
        N = cross(B, T);
    }

    frame.x = T;   // x = forward
    frame.y = N;   // y = left/right
    frame.z = B;   // z = up/down
}

Mat3f SL::FrameX(const Spline3& spline, float t)
{
    Vec3f v = Velocity(spline, t);

    Vec3f T = norm_safe(v);             // tangent
    Vec3f N = norm_safe(-cross_z(T));
    Vec3f B = cross(T, N);

    Mat3f frame;

    frame.x = T;   // x = forward
    frame.y = N;   // y = left
    frame.z = B;   // z = up

    return frame;
}

Mat3f SL::FrameX(const Spline3& spline, float t, Vec3f up)
{
    Vec3f v = Velocity(spline, t);

    Vec3f T = norm_safe(v);             // tangent
    Vec3f N = norm_safe(cross(up, T));  // force perpendicular to up
    Vec3f B = cross(T, N);

    Mat3f frame;
    frame.x = T;   // x = forward
    frame.y = N;   // y = left/right
    frame.z = B;   // z = up/down

    return frame;
}

float SL::LengthEstimate(const Spline3& s, float* error)
{
    // Our convex hull is p0, p1, p2, p3, so p0_p3 is our minimum possible length, and p0_p1 + p1_p2 + p2_p3 our maximum.
    float d03 = sqr(s.xb.x - s.xb.w) + sqr(s.yb.x - s.yb.w) + sqr(s.zb.x - s.zb.w);

    float d01 = sqr(s.xb.x - s.xb.y) + sqr(s.yb.x - s.yb.y) + sqr(s.zb.x - s.zb.y);
    float d12 = sqr(s.xb.y - s.xb.z) + sqr(s.yb.y - s.yb.z) + sqr(s.zb.y - s.zb.z);
    float d23 = sqr(s.xb.z - s.xb.w) + sqr(s.yb.z - s.yb.w) + sqr(s.zb.z - s.zb.w);

    float minLength = sqrtf(d03);
    float maxLength = sqrtf(d01) + sqrtf(d12) + sqrtf(d23);

    minLength *= 0.5f;
    maxLength *= 0.5f;

    *error = maxLength - minLength;
    return minLength + maxLength;
}

float SL::Length(const Spline3& s, float maxError)
{
    float error;
    float length = LengthEstimate(s, &error);

    if (error > maxError)
    {
        Spline3 s0;
        Spline3 s1;

        Split(s, &s0, &s1);

        return Length(s0, maxError) + Length(s1, maxError);
    }

    return length;
}

float SL::Length(const Spline3& s, float t0, float t1, float maxError)
{
    SL_ASSERT(t0 >= 0.0f && t0 <= 1.0f);
    SL_ASSERT(t1 >= 0.0f && t1 <= 1.0f);
    SL_ASSERT(t0 <= t1);

    Spline3 s0, s1;

    if (t0 == 0.0f)
    {
        if (t1 == 1.0f)
            return Length(s, maxError);

        Split(s, t1, &s0, &s1);
        return Length(s0, maxError);
    }
    else if (t1 >= 1.0f)
    {
        Split(s, t0, &s0, &s1);
        return Length(s1, maxError);
    }
    else
    {
        Split(s, t0, &s0, &s1);

        if (t1 == 1.0f)
            return Length(s1, maxError);

        Split(s1, (t1 - t0) / (1.0f - t0), &s0, &s1);
        return Length(s0, maxError);
    }
}

Bounds3 SL::FastBounds(const Spline3& spline)
{
    Vec2f bx = FastBounds(spline.xb);
    Vec2f by = FastBounds(spline.yb);
    Vec2f bz = FastBounds(spline.zb);

    return { { bx.x, by.x, bz.x }, { bx.y, by.y, bz.y } };
}

Bounds3 SL::ExactBounds(const Spline3& spline)
{
    Vec2f bx = ExactBounds(spline.xb);
    Vec2f by = ExactBounds(spline.yb);
    Vec2f bz = ExactBounds(spline.zb);

    return { { bx.x, by.x, bz.x }, { bx.y, by.y, bz.y } };
}


// Line conversion

int SL::SplinesToLinesLinear(int numSplines, const Spline3 splines[], tEmitLinesFunc3 emitLines, void* context, float step)
{
    Vec3f lineBuffer[kMaxEmitLines][2];
    int lineBufferCount = 0;
    int count = 0;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        lineBuffer[lineBufferCount][0] = Position0(spline);

        for (float t = step; t < 1.0f; t += step)
        {
            lineBuffer[lineBufferCount][1] = Position(spline, t);

            lineBufferCount++;
            if (lineBufferCount == kMaxEmitLines)
            {
                emitLines(lineBufferCount, lineBuffer, context);
                count += lineBufferCount;
                lineBufferCount = 0;
            }
            else
                lineBuffer[lineBufferCount][0] = lineBuffer[lineBufferCount - 1][1];
        }

        lineBuffer[lineBufferCount][1] = Position1(spline);

        lineBufferCount++;
        if (lineBufferCount == kMaxEmitLines)
        {
            emitLines(lineBufferCount, lineBuffer, context);
            count += lineBufferCount;
            lineBufferCount = 0;
        }
    }

    emitLines(lineBufferCount, lineBuffer, context);
    count += lineBufferCount;

    return count;
}

int SL::SplinesToLinesLinear(int numSplines, const Spline3 splines[], tEmitLinesExFunc3 emitLines, void* context, float step)
{
    Vec3f lineBuffer[kMaxEmitLines][2];
    float tvalBuffer[kMaxEmitLines][2];
    int lineBufferCount = 0;
    int count = 0;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        lineBuffer[lineBufferCount][0] = Position0(spline);
        tvalBuffer[lineBufferCount][0] = 0.0f;

        for (float t = step; t < 1.0f; t += step)
        {
            lineBuffer[lineBufferCount][1] = Position(spline, t);
            tvalBuffer[lineBufferCount][1] = t;

            int lastLineBufferCount = lineBufferCount++;

            if (lineBufferCount == kMaxEmitLines)
            {
                emitLines(lineBufferCount, lineBuffer, spline, i, numSplines, tvalBuffer, context);
                count += lineBufferCount;
                lineBufferCount = 0;
            }

            lineBuffer[lineBufferCount][0] = lineBuffer[lastLineBufferCount][1];
            tvalBuffer[lineBufferCount][0] = tvalBuffer[lastLineBufferCount][1];
        }

        lineBuffer[lineBufferCount][1] = Position1(spline);
        tvalBuffer[lineBufferCount][1] = 1.0f;

        lineBufferCount++;
        emitLines(lineBufferCount, lineBuffer, spline, i, numSplines, tvalBuffer, context);
        count += lineBufferCount;
        lineBufferCount = 0;
    }

    return count;
}

namespace
{
    float LineError(const Spline3& s)  // Calculate line approx error as max distance of hull points p1/p2 from p0_p1
    {
        Vec3f p0 = { s.xb.x, s.yb.x, s.zb.x };
        Vec3f p1 = { s.xb.y, s.yb.y, s.zb.y };
        Vec3f p2 = { s.xb.z, s.yb.z, s.zb.z };
        Vec3f p3 = { s.xb.w, s.yb.w, s.zb.w };

        Vec3f w = p3 - p0;
        float dww = dot(w, w);

        Vec3f v1 = p1 - p0;
        float dvw1 = dot(v1, w);
        Vec3f v2 = p2 - p0;
        float dvw2 = dot(v2, w);

        float error;

        if (dvw1 <= 0.0f)
            error = sqrlen(v1);
        else if (dvw1 >= dww)
            error = sqrlen(p1 - p3);
        else
            error = sqrlen((dvw1 / dww) * w - v1);

        if (dvw2 <= 0.0f)
            return vl_max(error, sqrlen(v2));
        else if (dvw2 >= dww)
            return vl_max(error, sqrlen(p2 - p3));
        else
            return vl_max(error, sqrlen((dvw2 / dww) * w - v2));
    }
}

int SL::SplinesToLinesAdaptive(int numSplines, const Spline3 splines[], tEmitLinesFunc3 emitLines, void* context, float tolerance)
{
    constexpr int kStackSize = 10; // gives us at least 1024 subdivs -- desirable to cap beyond this anyway.

    Spline3 stack [kStackSize];
    int      stackTop = 0;

    Vec3f lineBuffer[kMaxEmitLines][2];
    int   lineBufferCount = 0;

    int count = 0;
    tolerance *= tolerance;

    for (int i = 0; i < numSplines; i++)
    {
        stackTop = 0;
        stack [stackTop] = splines[i];

        while (stackTop >= 0)
        {
            SL_ASSERT_INDEX(stackTop, kStackSize);
            float error;

            if (stackTop < kStackSize - 1)
                error = LineError(stack[stackTop]);
            else
                error = 0.0f;

            if (error <= tolerance)
            {
                lineBuffer[lineBufferCount][0] = Position0(stack[stackTop]);
                lineBuffer[lineBufferCount][1] = Position1(stack[stackTop]);

                lineBufferCount++;
                if (lineBufferCount == kMaxEmitLines)
                {
                    emitLines(lineBufferCount, lineBuffer, context);
                    lineBufferCount = 0;
                }

                count++;
                stackTop--;
                continue;
            }

            Split(stack[stackTop], stack + stackTop + 1, stack + stackTop + 0);
            stackTop++;
        }
    }

    if (lineBufferCount > 0)
        emitLines(lineBufferCount, lineBuffer, context);

    return count;
}

int SL::SplinesToLinesAdaptive(int numSplines, const Spline3 splines[], tEmitLinesExFunc3 emitLines, void* context, float tolerance)
{
    constexpr int kStackSize = 10; // gives us at least 1024 subdivs -- desirable to cap beyond this anyway.

    Spline3 stack [kStackSize];
    Vec2f    stackT[kStackSize];   // for start/end t values
    int      stackTop = 0;

    Vec3f lineBuffer[kMaxEmitLines][2];
    float tvalBuffer[kMaxEmitLines][2];
    int   lineBufferCount = 0;

    int count = 0;
    tolerance *= tolerance;

    for (int i = 0; i < numSplines; i++)
    {
        stackTop = 0;
        stack [stackTop] = splines[i];
        stackT[stackTop] = { 0.0f, 1.0f };

        while (stackTop >= 0)
        {
            SL_ASSERT_INDEX(stackTop, kStackSize);
            float error;

            if (stackTop < kStackSize - 1)
                error = LineError(stack[stackTop]);
            else
                error = 0.0f;

            if (error <= tolerance)
            {
                lineBuffer[lineBufferCount][0] = Position0(stack[stackTop]);
                lineBuffer[lineBufferCount][1] = Position1(stack[stackTop]);

                tvalBuffer[lineBufferCount][0] = stackT[stackTop].x;
                tvalBuffer[lineBufferCount][1] = stackT[stackTop].y;

                lineBufferCount++;
                if (lineBufferCount == kMaxEmitLines)
                {
                    emitLines(lineBufferCount, lineBuffer, splines[i], i, numSplines, tvalBuffer, context);
                    lineBufferCount = 0;
                }

                count++;
                stackTop--;
                continue;
            }

            Split(stack[stackTop], stack + stackTop + 1, stack + stackTop + 0);

            float halfT = 0.5f * (stackT[stackTop].x + stackT[stackTop].y);
            stackT[stackTop + 1].x = stackT[stackTop + 0].x;
            stackT[stackTop + 1].y = halfT;
            stackT[stackTop + 0].x = halfT;
            stackTop++;
        }

        if (lineBufferCount > 0)
        {
            emitLines(lineBufferCount, lineBuffer, splines[i], i, numSplines, tvalBuffer, context);
            lineBufferCount = 0;
        }
    }

    return count;
}

// Subdivision

void SL::Split(const Spline3& spline, Spline3* spline0, Spline3* spline1)
{
    Split(spline.xb, &spline0->xb, &spline1->xb);
    Split(spline.yb, &spline0->yb, &spline1->yb);
    Split(spline.zb, &spline0->zb, &spline1->zb);
}

void SL::Split(const Spline3& spline, float t, Spline3* spline0, Spline3* spline1)
{
    Split(spline.xb, t, &spline0->xb, &spline1->xb);
    Split(spline.yb, t, &spline0->yb, &spline1->yb);
    Split(spline.zb, t, &spline0->zb, &spline1->zb);
}

bool SL::Join(const Spline3& s0, const Spline3& s1, Spline3* splineOut)
{
    return
       Join(s0.xb, s1.xb, &splineOut->xb)
    && Join(s0.yb, s1.yb, &splineOut->yb)
    && Join(s0.zb, s1.zb, &splineOut->zb);
}

void SL::Split(vector<Spline3>* splinesIn)
{
    vector<Spline3> splines;

    for (const Spline3& s : *splinesIn)
    {
        Spline3 s0, s1;

        Split(s, &s0, &s1);
        splines.push_back(s0);
        splines.push_back(s1);
    }

    splinesIn->swap(splines);
}

void SL::Split(vector<Spline3>* splinesIn, int n)
{
    vector<Spline3> splines;

    for (const Spline3& s : *splinesIn)
    {
        Spline3 ss(s);
        Spline3 s0, s1;

        for (int i = n; i > 1; i--)
        {
            Split(ss, 1.0f / i, &s0, &ss);
            splines.push_back(s0);
        }
        splines.push_back(ss);
    }

    splinesIn->swap(splines);
}

void SL::Join(vector<Spline3>* splinesIn)
{
    vector<Spline3> splines;
    const Spline3* prevS = 0;

    for (const Spline3& s : *splinesIn)
    {
        if (!prevS)
        {
            prevS = &s;
            continue;
        }

        Spline3 sj;
        if (Join(*prevS, s, &sj))
            splines.push_back(sj);
        else
        {
            splines.push_back(*prevS);
            splines.push_back(s);
        }

        prevS = 0;
    }

    if (prevS)
        splines.push_back(*prevS);

    splinesIn->swap(splines);
}

Spline3 SL::Trim(const Spline3& spline, float t0, float t1)
{
    return { Trim(spline.xb, t0, t1), Trim(spline.yb, t0, t1), Trim(spline.zb, t0, t1) };
}

namespace
{
    void SubdivideForLength(const Spline3& s, vector<Spline3>* splines, float tolerance)
    {
        float error;
        float length = LengthEstimate(s, &error);

        if (error <= tolerance * length)
            splines->push_back(s);
        else
        {
            Spline3 s1, s2;
            Split(s, &s1, &s2);

            SubdivideForLength(s1, splines, tolerance);
            SubdivideForLength(s2, splines, tolerance);
        }
    }
}

void SL::SubdivideForLength(vector<Spline3>* splinesIn, float tolerance)
{
    vector<Spline3> splines;

    for (const Spline3& s : *splinesIn)
        ::SubdivideForLength(s, &splines, tolerance);

    splinesIn->swap(splines);
}

namespace
{
    float ArcError(const Spline3& s, float* tSplit)
    {
        Vec2f ex = ArcError2(s.xb);
        Vec2f ey = ArcError2(s.yb);
        Vec2f ez = ArcError2(s.zb);

        float e0 = ex.x + ey.x + ez.x;
        float e1 = ex.y + ey.y + ez.y;
        float es2 = e0 + e1;

        float f = (es2 < 1e-6f) ? 0.5f : sqrtf(e0 / es2);

        *tSplit = (1.0f / 3.0f) * (1.0f + f);

        return sqrtf(es2);
    }

    void SubdivideForT(const Spline3& s, vector<Spline3>* splines, float tolerance)
    {
        float splitT;
        float err = ArcError(s, &splitT);

        if (err <= tolerance)
            splines->push_back(s);
        else
        {
            Spline3 s1, s2;
            Split(s, splitT, &s1, &s2);

            SubdivideForT(s1, splines, tolerance);
            SubdivideForT(s2, splines, tolerance);
        }
    }
}

void SL::SubdivideForT(vector<Spline3>* splinesIn, float tolerance)
{
    vector<Spline3> splines;

    for (const Spline3& s : *splinesIn)
        ::SubdivideForT(s, &splines, tolerance);

    splinesIn->swap(splines);
}

void SL::SubdivideForLengthRatio(vector<Vec3f>& positions, float maxRatio)
{
    for (int i = 1; i < size_i(positions); )
    {
        float l1 = len(positions[i] - positions[i - 1]);

        if (i > 1)
        {
            float l0 = len(positions[i - 1] - positions[i - 2]);

            if (l1 > l0 * maxRatio)
            {
                positions.insert(positions.begin() + i, 0.5f * (positions[i] + positions[i - 1]));
                continue;
            }
        }

        if (i < size_i(positions) - 1)
        {
            float l2 = len(positions[i + 1] - positions[i]);

            if (l1 > l2 * maxRatio)
            {
                positions.insert(positions.begin() + i, 0.5f * (positions[i] + positions[i - 1]));
                continue;
            }
        }

        i++;
    }
}

namespace
{
    inline float ClosestPoint(Vec3f p, Vec3f p0, Vec3f p1)
    {
        Vec3f w = p1 - p0;
        Vec3f v =  p - p0;

        float dvw = dot(v, w);

        if (dvw <= 0.0f)
            return 0.0f;

        float dww = dot(w, w);

        if (dvw >= dww)
            return 1.0f;

        return dvw / dww;
    }

    void FindClosestPointNewtonRaphson(const Spline3& spline, Vec3f p, float sIn, int maxIterations, float* tOut, float* dOut)
    {
        SL_ASSERT(sIn >= 0.0f && sIn <= 1.0f);
        const float maxS = 1.0f - 1e-6f;

        float sk = sIn;
        float dk = len(Position(spline, sk) - p);

        constexpr float width = 1e-3f;

        float maxJump  = 0.5f;   // avoid jumping too far, leads to oscillation

        for (int i = 0; i < maxIterations; i++)
        {
            float ss = Clamp(sk, width, 1.0f - width); // so can interpolate points for Newtons method

            float d1 = len(Position(spline, ss - width) - p);
            float d2 = len(Position(spline, ss        ) - p);
            float d3 = len(Position(spline, ss + width) - p);

            float g1 = (d2 - d1) / width;
            float g2 = (d3 - d2) / width;

            float grad = (d3 - d1) / (2.0f * width);
            float curv = (g2 - g1) / width;

            float sn;   // next candidate

            if (curv > 0.0f)    // if d' is heading towards a minima, apply NR for d'
                sn = ss - grad / curv;
            else if (grad != 0.0f)
                sn = ss - d2 / grad; // otherwise, apply for D.
            else
                sn = sk;

            sn = Clamp(sn, sk - maxJump, sk + maxJump);   // avoid large steps, often unstable.

            // only update our estimate if the new value is in range and closer.
            if (sn >= 0.0f && sn < maxS)
            {
                float dn = len(Position(spline, sn) - p);

                if (dn < dk)    // only update sk if d actually gets smaller
                {
                    sk = sn;
                    dk = dn;
                }
            }

            maxJump *= 0.5f;    // reduce on a schedule -- helps binary search back towards a jump that is valid.
        }

        (*tOut) = sk;
        (*dOut) = dk;
    }
}

float SL::FindClosestPoint(Vec3f p, const Spline3& spline)
{
    // Approximate s from straight line between the start and end.
    float s = ClosestPoint(p, Position0(spline), Position1(spline));

    // Use this as starting point for Newton-Raphson solve.
    float d;
    FindClosestPointNewtonRaphson(spline, p, s, 8, &s, &d);

    return s;
}

float SL::FindClosestPoint(Vec3f p, int numSplines, const Spline3 splines[], int* index)
{
    vector<SubSpline3> nearbyInfo;

    FindNearbySplines(p, numSplines, splines, &nearbyInfo);
    return FindClosestPoint(p, numSplines, splines, nearbyInfo, index);
}


namespace
{
    void FindMinMaxDistance2s(Vec3f p, const Bounds3& bbox, float* minD2, float* maxD2)
    {
        const Vec3f& p0 = bbox.mMin;
        const Vec3f& p1 = bbox.mMax;

        // Find the nearest point to p inside the bbox
        // This can be a bbox vertex, a point on an edge or face, or p itself if it's inside the box
        float minX = Clamp(p.x, p0.x, p1.x);
        float minY = Clamp(p.y, p0.y, p1.y);
        float minZ = Clamp(p.z, p0.z, p1.z);

        // Find the farthest point from p inside the bbox
        // This is always a bbox vertex.
        Vec3f d0(abs(p - p0));
        Vec3f d1(abs(p - p1));

        float maxX = d0.x > d1.x ? p0.x : p1.x; // select the coordinate we're farthest from
        float maxY = d0.y > d1.y ? p0.y : p1.y;
        float maxZ = d0.z > d1.z ? p0.z : p1.z;

        // return len2
        *minD2 = sqr(p.x - minX) + sqr(p.y - minY) + sqr(p.z - minZ);
        *maxD2 = sqr(p.x - maxX) + sqr(p.y - maxY) + sqr(p.z - maxZ);
    }

    void FindMinMaxDistance2s(Vec3f p, const Spline3& spline, float* minD2, float* maxD2)
    {
        Bounds3 bbox = FastBounds(spline);

        FindMinMaxDistance2s(p, bbox, minD2, maxD2);
    }

    void Split(const SubSpline3& s, SubSpline3* s0, SubSpline3* s1)
    {
        ::Split(s.mSpline.xb, &s0->mSpline.xb, &s1->mSpline.xb);
        ::Split(s.mSpline.yb, &s0->mSpline.yb, &s1->mSpline.yb);
        ::Split(s.mSpline.zb, &s0->mSpline.zb, &s1->mSpline.zb);

        s0->mParent = s.mParent;
        s1->mParent = s.mParent;
    }
}

int SL::FindNearbySplines(Vec3f p, int numSplines, const Spline3 splines[], vector<SubSpline3>* results, float* smallestFarOut, int numIter)
{
    vector<SubSpline3>& nearSplines = *results;

    nearSplines.clear();

    float smallestFar  = FLT_MAX;
    float smallestNear = FLT_MAX;

    // Find initial list

    int maxSize = 0;

    for (int i = 0; i < numSplines; i++)
    {
        float near;
        float far;
        FindMinMaxDistance2s(p, splines[i], &near, &far);

        if (near < smallestFar)
        {
            // we at least overlap the current set.
            if (near < smallestNear)
                smallestNear = near;

            if (far < smallestFar)
            {
                // we bring in the 'best' far distance
                smallestFar = far;

                // compact list to reject any segments that now cannot be closest.
                int dj = 0;

                for (int j = 0, nj = size_i(nearSplines); j < nj; j++)
                    if (nearSplines[j].mD2 <= smallestFar)
                    {
                        if (dj < j)
                            nearSplines[dj] = nearSplines[j];

                        dj++;
                    }

                nearSplines.resize(dj);
            }

            SubSpline3 ss = { splines[i], i, near };
            nearSplines.push_back(ss);

            if (maxSize < size_i(nearSplines))
                maxSize = size_i(nearSplines);
        }
    }

    // Subdivide + refine

    int numNearSplines = size_i(nearSplines);

    for (int i = 0; i < numIter; i++)
    {
        int numNearSplines2 = numNearSplines * 2;

        nearSplines.resize(numNearSplines2);

        for (int i = numNearSplines - 1; i >= 0; i--)
            ::Split(nearSplines[i], &nearSplines[2 * i], &nearSplines[2 * i + 1]);

        smallestNear = FLT_MAX; // this may actually increase on subdivision.

        for (int i = 0; i < numNearSplines2; i++)
        {
            float near;
            float far;
            FindMinMaxDistance2s(p, nearSplines[i].mSpline, &near, &far);

            if (far < smallestFar)
                smallestFar = far;
            if (near < smallestNear)
                smallestNear = near;

            nearSplines[i].mD2 = near;
        }

        int di = 0;
        for (int i = 0; i < numNearSplines2; i++)
            if (nearSplines[i].mD2 <= smallestFar)
            {
                if (di < i)
                    nearSplines[di] = nearSplines[i];

                di++;
            }

        nearSplines.resize(di);
        numNearSplines = di;
    }

    SL_ASSERT(numNearSplines > 0 || numSplines == 0);

    if (smallestFarOut)
        *smallestFarOut = smallestFar;

    return numNearSplines;
}

float SL::FindClosestPoint(Vec3f p, int numSplines, const Spline3 splines[], const vector<SubSpline3>& nearbySplines, int* index)
{
    int prevParent = -1;
    float minD = FLT_MAX;
    float minT = 0.0f;

    *index = -1;

    for (const SubSpline3& subSpline : nearbySplines)
    {
        if (subSpline.mParent != prevParent)
        {
            SL_ASSERT(subSpline.mParent >= 0 && subSpline.mParent < numSplines);
            (void) numSplines;

            const Spline3& spline = splines[subSpline.mParent];

            float t = ClosestPoint(p, Position0(spline), Position1(spline));
            float d;

            FindClosestPointNewtonRaphson(spline, p, t, 8, &t, &d);

            if (minD > d)
            {
                minD = d;
                *index = subSpline.mParent;
                minT = t;
            }
        }
    }

    return minT;
}

namespace
{
    struct SubSplineT3
    {
        Spline3 mSpline;
        float    mT0;
        float    mT1;
    };

    inline void Split(const SubSplineT3& s, SubSplineT3* s0, SubSplineT3* s1)
    {
        ::Split(s.mSpline.xb, &s0->mSpline.xb, &s1->mSpline.xb);
        ::Split(s.mSpline.yb, &s0->mSpline.yb, &s1->mSpline.yb);
        ::Split(s.mSpline.zb, &s0->mSpline.zb, &s1->mSpline.zb);

        float midT = 0.5f * (s.mT0 + s.mT1);

        s0->mT0     = s.mT0;
        s0->mT1     = midT;

        s1->mT0     = midT;
        s1->mT1     = s.mT1;
    }

    int FindSubSplineIntersections
    (
        const SubSplineT3& spline0,
        const SubSplineT3& spline1,
        int   dest,
        int   maxDest,
        float results[][2],
        float tolerance
    )
    {
        SL_ASSERT(dest < maxDest);

        Bounds3 bbox0 = ExactBounds(spline0.mSpline);
        Bounds3 bbox1 = ExactBounds(spline1.mSpline);

        if (!Intersects(bbox0, bbox1))
            return dest;

        if (Larger(bbox0, tolerance))
        {
            SubSplineT3 spline00, spline01;
            Split(spline0, &spline00, &spline01);

            if (Larger(bbox1, tolerance))
            {
                SubSplineT3 spline10, spline11;
                Split(spline1, &spline10, &spline11);

                dest     = FindSubSplineIntersections(spline00, spline10, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline01, spline10, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline00, spline11, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline01, spline11, dest, maxDest, results, tolerance);
            }
            else
            {
                dest     = FindSubSplineIntersections(spline00, spline1, dest, maxDest, results, tolerance);
                if (dest < maxDest)
                    dest = FindSubSplineIntersections(spline01, spline1, dest, maxDest, results, tolerance);
            }

            return dest;
        }

        if (Larger(bbox1, tolerance))
        {
            SubSplineT3 spline10, spline11;
            Split(spline1, &spline10, &spline11);

            dest     = FindSubSplineIntersections(spline0, spline10, dest, maxDest, results, tolerance);
            if (dest < maxDest)
                dest = FindSubSplineIntersections(spline0, spline11, dest, maxDest, results, tolerance);

            return dest;
        }

        float t0 = 0.5f * (spline0.mT0 + spline0.mT1);
        float t1 = 0.5f * (spline1.mT0 + spline1.mT1);

        // debounce
        for (int i = 0; i < dest; i++)
            if (fabsf(results[i][0] - t0) < 1e-2f
             && fabsf(results[i][1] - t1) < 1e-2f)
            {
                return dest;
            }

        results[dest][0] = t0;
        results[dest][1] = t1;
        return dest + 1;
    }
}

int SL::FindSplineIntersections(const Spline3& spline0, const Spline3& spline1, int maxResults, float results[][2], float tolerance)
{
    if (maxResults <= 0)
        return 0;

    SubSplineT3 subSpline0 = { spline0, 0.0f, 1.0f };
    SubSplineT3 subSpline1 = { spline1, 0.0f, 1.0f };

    return FindSubSplineIntersections(subSpline0, subSpline1, 0, maxResults, results, tolerance);
}

namespace
{
    int FindSplineIntersections
    (
        int   is0, int numSplines0, const Spline3 splines0[],
        int   is1, int numSplines1, const Spline3 splines1[],
        int   maxResults,
        int   resultsI[][2],
        float resultsT[][2],
        float tolerance
    )
    {
        if (maxResults <= 0 || numSplines0 == 0 || numSplines1 == 0)
            return 0;

        int count = 0;

        // once the lists are small enough, brute-force the remainder, as recalculating the bounds is not free
        if (numSplines0 >= 4 || numSplines1 >= 4)
        {
            // Terminate if the lists don't overlap
            Bounds3 b0 = FastBounds(splines0[0]);
            Bounds3 b1 = FastBounds(splines1[0]);

            for (int i = 1; i < numSplines0; i++)
                Add(b0, FastBounds(splines0[i]));
            for (int i = 1; i < numSplines1; i++)
                Add(b1, FastBounds(splines1[i]));

            if (!Intersects(b0, b1))
                return 0;

            // Divide each spline list into two, and recurse
            int n00 = numSplines0 / 2;
            int n01 = numSplines0 - n00;
            int n10 = numSplines1 / 2;
            int n11 = numSplines1 - n10;

            count += FindSplineIntersections(is0 +   0, n00, splines0 +   0, is1 +   0, n10, splines1 +   0, maxResults - count, resultsI + count, resultsT + count, tolerance);
            count += FindSplineIntersections(is0 +   0, n00, splines0 +   0, is1 + n10, n11, splines1 + n10, maxResults - count, resultsI + count, resultsT + count, tolerance);
            count += FindSplineIntersections(is0 + n00, n01, splines0 + n00, is1 +   0, n10, splines1 +   0, maxResults - count, resultsI + count, resultsT + count, tolerance);
            count += FindSplineIntersections(is0 + n00, n01, splines0 + n00, is1 + n10, n11, splines1 + n10, maxResults - count, resultsI + count, resultsT + count, tolerance);

            return count;
        }

        SubSplineT3 st0, st1;

        st0.mT0 = 0.0f;
        st0.mT1 = 1.0f;
        st1.mT0 = 0.0f;
        st1.mT1 = 1.0f;

        for (int i0 = 0; i0 < numSplines0; i0++)
            for (int i1 = 0; i1 < numSplines1; i1++)
            {
                st0.mSpline = splines0[i0];
                st1.mSpline = splines1[i1];

                int numIntersections = FindSubSplineIntersections(st0, st1, 0, maxResults - count, resultsT + count, tolerance);

                for (int k = 0; k < numIntersections; k++)
                {
                    resultsI[k + count][0] = is0 + i0;
                    resultsI[k + count][1] = is1 + i1;
                }

                count += numIntersections;
                SL_ASSERT(count <= maxResults);

                if (count == maxResults)
                    return count;
            }

        return count;
    }
}

int SL::FindSplineIntersections
(
    int   numSplines0, const Spline3 splines0[],
    int   numSplines1, const Spline3 splines1[],
    int   maxResults,
    int   resultsI[][2],
    float resultsT[][2],
    float tolerance
)
{
    return ::FindSplineIntersections(0, numSplines0, splines0, 0, numSplines1, splines1, maxResults, resultsI, resultsT, tolerance);
}

int SL::FindSplineIntersections
(
    int   numSplines, const Spline3 splines[],
    int   maxResults,
    int   resultsI[][2],
    float resultsT[][2],
    float tolerance
)
// Find self-intersections
{
    if (maxResults <= 0 || numSplines == 0)
        return 0;

    int count = 0;

    if (numSplines >= 8)
    {
        const int n0 = numSplines / 2;
        const int n1 = numSplines - n0;

        // Find intersections between the two halves
        count += ::FindSplineIntersections(0, n0, splines, 0, n1, splines + n0, maxResults - count, resultsI + count, resultsT + count, tolerance);

        for (int i = 0; i < count; i++)
        {
            // ignore spurious intersections between endpoint of first and start of second
            if (resultsI[i][1] == 0 && resultsI[i][0] == n0 - 1)
                if (resultsT[i][0] > 0.95f && resultsT[i][1] < 0.05f)
                {
                    resultsT[i][0] = resultsT[count - 1][0];
                    resultsT[i][1] = resultsT[count - 1][1];
                    resultsI[i][0] = resultsI[count - 1][0];
                    resultsI[i][1] = resultsI[count - 1][1];

                    i--;
                    count--;
                    continue;
                }

            resultsI[i][1] += n0;           // adjust second half indices
        }

        // Find self-intersections in the first half
        count += FindSplineIntersections(n0, splines +  0, maxResults - count, resultsI + count, resultsT + count, tolerance);

        // Find self-intersections in the second half
        int prevCount = count;
        count += FindSplineIntersections(n1, splines + n0, maxResults - count, resultsI + count, resultsT + count, tolerance);

        for (int i = prevCount; i < count; i++)
        {
            resultsI[i][0] += n0;   // adjust second half indices
            resultsI[i][1] += n0;
        }

        return count;
    }

    SubSplineT3 st0, st1;

    st0.mT0 = 0.0f;
    st0.mT1 = 1.0f;
    st1.mT0 = 0.0f;
    st1.mT1 = 1.0f;

    for (int i0 = 0; i0 < numSplines; i0++)
        for (int i1 = i0 + 1; i1 < numSplines; i1++)
        {
            st0.mSpline = splines[i0];
            st1.mSpline = splines[i1];

            int numIntersections = FindSubSplineIntersections(st0, st1, 0, maxResults, resultsT, tolerance);

            if (i1 == i0 + 1)   // ignore spurious intersections between endpoint of i0 and start of i1
            {
                for (int k = 0; k < numIntersections; k++)
                {
                    if (resultsT[k][0] > 0.95f && resultsT[k][1] < 0.05f)
                    {
                        resultsT[k][0] = resultsT[numIntersections - 1][0];
                        resultsT[k][1] = resultsT[numIntersections - 1][1];

                        k--;
                        numIntersections--;
                    }
                }
            }

            for (int k = 0; k < numIntersections; k++)
            {
                resultsI[k][0] = i0;
                resultsI[k][1] = i1;
            }

            count += numIntersections;
            maxResults -= numIntersections;

            if (maxResults <= 0)
                return count;

            resultsT += numIntersections;
            resultsI += numIntersections;
        }

    return count;
}

float SL::AdvanceAgent(const Spline3& spline, float t, float ds)
{
    Vec3f v  = Velocity(spline, t);
    float v2 = sqrlen(v);
    float dt = ds;

    if (v2 > 0.01f)
        dt *= InvSqrtFast(v2);
    else
        dt *= 10.0f;

    return t + dt;
}

bool SL::AdvanceAgent(int* index, float* t, int numSplines, const Spline3 splines[], float ds)
{
    *t = AdvanceAgent(splines[*index], *t, ds);

    return ::AdvanceAgent(index, t, numSplines);
}

Spline3 SL::Reverse(const Spline3& spline)
{
    return
    {
        Vec4f(spline.xb.w, spline.xb.z, spline.xb.y, spline.xb.x),
        Vec4f(spline.yb.w, spline.yb.z, spline.yb.y, spline.yb.x),
        Vec4f(spline.zb.w, spline.zb.z, spline.zb.y, spline.zb.x)
    };
}

void SL::Reverse(vector<Spline3>* splines)
{
    int n = int(splines->size());
    int h = n / 2;

    for (int i = 0; i < h; i++)
    {
        Spline3& s0 = (*splines)[i];
        Spline3& s1 = (*splines)[n - i - 1];

        Spline3 sr0 = Reverse(s1);
        Spline3 sr1 = Reverse(s0);

        s0 = sr0;
        s1 = sr1;
    }

    if (2 * h < n)
        (*splines)[h] = Reverse((*splines)[h]);
}

Spline3 SL::Offset(const Spline3& spline, float offset)
{
    float sx0 = spline.xb.y - spline.xb.x;
    float sy0 = spline.yb.y - spline.yb.x;
    float sd0 = InvSqrtFast(sx0 * sx0 + sy0 * sy0) * offset;

    float sx1 = spline.xb.z - spline.xb.x;
    float sy1 = spline.yb.z - spline.yb.x;
    float sd1 = InvSqrtFast(sx1 * sx1 + sy1 * sy1) * offset;

    float sx2 = spline.xb.w - spline.xb.y;
    float sy2 = spline.yb.w - spline.yb.y;
    float sd2 = InvSqrtFast(sx2 * sx2 + sy2 * sy2) * offset;

    float sx3 = spline.xb.w - spline.xb.z;
    float sy3 = spline.yb.w - spline.yb.z;
    float sd3 = InvSqrtFast(sx3 * sx3 + sy3 * sy3) * offset;

    return
    {
        spline.xb + Vec4f(sy0 * sd0, sy1 * sd1, sy2 * sd2, sy3 * sd3),
        spline.yb - Vec4f(sx0 * sd0, sx1 * sd1, sx2 * sd2, sx3 * sd3),
        spline.zb
    };
}

void SL::Offset(vector<Spline3>* splines, float offset)
{
    for (Spline3& s : *splines)
        s = Offset(s, offset);
}
