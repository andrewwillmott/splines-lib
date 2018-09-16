//
//  File:       Splines.cpp
//
//  Function:   Various cubic spline utilities
//
//  Copyright:  Andrew Willmott 2018
//

#include "Splines.h"

#include <float.h>
#include <math.h>

using namespace SplineLib;

namespace
{
    // Mini VL
    inline float sqr(float x) { return x * x; }
    inline float lerp(float a, float b, float t) { return (1.0f - t) * a + t * b; }

    inline Vec2f operator+(Vec2f a, Vec2f b) { return { a.x + b.x, a.y + b.y }; }
    inline Vec2f operator-(Vec2f a, Vec2f b) { return { a.x - b.x, a.y - b.y }; }
    inline Vec2f operator*(float s, Vec2f a) { return { s * a.x, s * a.y }; }
    inline float dot(Vec2f a, Vec2f b)       { return a.x * b.x + a.y * b.y; }
    inline float len      (Vec2f v)          { return sqrtf(v.x * v.x + v.y * v.y); }
    inline float sqrlen   (Vec2f v)          { return v.x * v.x + v.y * v.y; }
    inline Vec2f norm_safe(Vec2f v)          { return (1.0f / (len(v) + 1e-8f)) * v; }
    inline Vec2f cross    (Vec2f v)          { return { -v.y, v.x }; }
    inline Vec2f abs      (Vec2f v)          { return { fabsf(v.x), fabsf(v.y) }; }

    inline Vec3f operator-(Vec3f v)          { return { -v.x, -v.y, -v.z }; }
    inline Vec3f operator+(Vec3f a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z}; }
    inline Vec3f operator-(Vec3f a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z}; }
    inline Vec3f operator*(float s, Vec3f a) { return { s   * a.x, s   * a.y, s   * a.z}; }
    inline float dot      (Vec3f a, Vec3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    inline float len      (Vec3f v)          { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); }
    inline float sqrlen   (Vec3f v)          { return v.x * v.x + v.y * v.y + v.z * v.z; }
    inline Vec3f norm     (Vec3f v)          { return (1.0f / len(v)) * v; }
    inline Vec3f norm_safe(Vec3f v)          { return (1.0f / (len(v) + 1e-8f)) * v; }
    inline Vec3f abs      (Vec3f v)          { return { fabsf(v.x), fabsf(v.y), fabsf(v.z) }; }
    inline Vec3f cross    (Vec3f a, Vec3f b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
    inline Vec3f cross_z  (Vec3f v)          { return { v.y, -v.x, 0.0f }; }

    inline Vec4f operator-(Vec4f v)          { return { -v.x, -v.y, -v.z, -v.w }; }
    inline Vec4f operator+(Vec4f a, Vec4f b) { return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w}; }
    inline Vec4f operator-(Vec4f a, Vec4f b) { return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w}; }
    inline Vec4f operator*(float s, Vec4f a) { return { s   * a.x, s   * a.y, s   * a.z, s   * a.w}; }
    inline float dot   (Vec4f a, Vec4f b)    { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

    inline float Max(float a, float b)       { return b < a ? a : b; }
    inline float Min(float a, float b)       { return a < b ? a : b; }

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
    
    inline bool Larger(const Bounds2f& bb, float t) { Vec2f d = bb.max - bb.min; return d.x > t || d.y > t; }
    inline bool Larger(const Bounds3f& bb, float t) { Vec3f d = bb.max - bb.min; return d.x > t || d.y > t || d.z < t; }
    inline bool Intersects(const Bounds2f& a, const Bounds2f& b)
    {
        return a.max.x >= b.min.x && a.min.x <= b.max.x
            && a.max.y >= b.min.y && a.min.y <= b.max.y;
    }
    inline bool Intersects(const Bounds3f& a, const Bounds3f& b)
    {
        return a.max.x >= b.min.x && a.min.x <= b.max.x
            && a.max.y >= b.min.y && a.min.y <= b.max.y
            && a.max.z >= b.min.z && a.min.z <= b.max.z;
    }
    inline void Add(Bounds2f& a, const Bounds2f& b)
    {
        if (a.min.x > b.min.x) a.min.x = b.min.x; else if (a.max.x < b.max.x) a.max.x = b.max.x;
        if (a.min.y > b.min.y) a.min.y = b.min.y; else if (a.max.y < b.max.y) a.max.y = b.max.y;
    }
    inline void Add(Bounds3f& a, const Bounds3f& b)
    {
        if (a.min.x > b.min.x) a.min.x = b.min.x; else if (a.max.x < b.max.x) a.max.x = b.max.x;
        if (a.min.y > b.min.y) a.min.y = b.min.y; else if (a.max.y < b.max.y) a.max.y = b.max.y;
        if (a.min.z > b.min.z) a.min.z = b.min.z; else if (a.max.z < b.max.z) a.max.z = b.max.z;
    }

    template<class T> inline int size_i(const T& container) { return int(container.size()); }
}

namespace
{
    // Utilities
    
    inline Vec4f BezierWeights(float t)
    /// Returns Bezier basis weights for 't'
    {
        float s  = 1.0f - t;
        
        float t2 = t * t;
        float t3 = t2 * t;
        
        float s2 = s * s;
        float s3 = s2 * s;
        
        return Vec4f(s3, 3.0f * s2 * t, 3.0f * s * t2, t3);
    }
    
    inline Vec4f BezierWeights(const Vec4f& t)
    /// Vector version useful for derivatives
    {
        return Vec4f
        (
            t.x - 3.0f * t.y + 3.0f * t.z -        t.w,
                  3.0f * t.y - 6.0f * t.z + 3.0f * t.w,
                               3.0f * t.z - 3.0f * t.w,
                                                   t.w
        );
    }

    inline Vec4f CubicCoeffs(const Vec4f& b)
    /// Returns cubic coefficients for the given Bezier weights
    {
        return Vec4f
        (
                    b.x                                ,
            -3.0f * b.x + 3.0f * b.y                   ,
             3.0f * b.x - 6.0f * b.y + 3.0f * b.z      ,
                   -b.x + 3.0f * b.y - 3.0f * b.z + b.w
        );
    }

    inline Vec2f HullBounds(const Vec4f& s)
    /// Returns bounds of the convex hull
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

        return Vec2f
        (
            Min(b01.x, b23.x),
            Max(b01.y, b23.y)
        );
    }

    Vec2f ExactBounds(const Vec4f& spline)
    /// Returns accurate bounds taking extrema into account.
    {
        Vec2f bounds;

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
    inline void Split(float t, const Vec4f& spline, Vec4f* spline0, Vec4f* spline1)
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

        *spline0 = Vec4f(sx, q0, r0, s0);
        *spline1 = Vec4f(s0, r1, q2, sw);
    }

    // Optimised for t=0.5
    inline void Split(const Vec4f& spline, Vec4f* spline0, Vec4f* spline1)
    {
        float q0 = (spline.x + spline.y) * 0.5f;    // x + y / 2
        float q1 = (spline.y + spline.z) * 0.5f;    // y + z / 2
        float q2 = (spline.z + spline.w) * 0.5f;    // z + w / 2

        float r0 = (q0 + q1) * 0.5f;    // x + 2y + z / 4
        float r1 = (q1 + q2) * 0.5f;    // y + 2z + w / 4

        float s0 = (r0 + r1) * 0.5f;    // q0 + 2q1 + q2 / 4 = x+y + 2(y+z) + z+w / 8 = x + 3y + 3z + w

        float sx = spline.x;    // support aliasing
        float sw = spline.w;

        *spline0 = Vec4f(sx, q0, r0, s0);
        *spline1 = Vec4f(s0, r1, q2, sw);
    }

    bool Join(const Vec4f& s0, const Vec4f& s1, Vec4f* sOut)
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
            *sOut = Vec4f(x0, y0, z1, w1);   // use most stable terms
            return true;
        }

        return false;
    }

    inline Vec2f ArcError2(Vec4f s)
    /// Returns squared displacement from linear (b0_b3) for hull points b1/b2
    {
        float w = s.w - s.x;

        float ty = s.x + w * 1.0f / 3.0f - s.y;
        float tz = s.x + w * 2.0f / 3.0f - s.z;
        float d2 = 1.0f / (sqr(w) + 1.0f);

        return Vec2f(sqr(ty) * d2, sqr(tz) * d2);
    }

    bool AdvanceAgent(int* indexInOut, float* tInOut, int numSplines)
    /// Update index for t if necessary, but don't run off array
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
// 2D
////////////////////////////////////////////////////////////////////////////////


cSpline2 SplineLib::BezierSpline(const Vec2f& p0, const Vec2f& p1, const Vec2f& p2, const Vec2f& p3)
{
    return cSpline2
    {
        Vec4f(p0.x, p1.x, p2.x, p3.x),
        Vec4f(p0.y, p1.y, p2.y, p3.y),
    };
}

cSpline2 SplineLib::HermiteSpline(const Vec2f& p0, const Vec2f& p1, const Vec2f& v0, const Vec2f& v1)
{
    Vec2f pb1 = p0 + (1.0f / 3.0f) * v0;
    Vec2f pb2 = p1 - (1.0f / 3.0f) * v1;

    return BezierSpline(p0, pb1, pb2, p1);
}

cSpline2 SplineLib::CatmullRomSpline(const Vec2f& p0, const Vec2f& p1, const Vec2f& p2, const Vec2f& p3)
{
    Vec2f pb1 = p1 + (1.0f / 6.0f) * (p2 - p0);
    Vec2f pb2 = p2 - (1.0f / 6.0f) * (p3 - p1);

    return BezierSpline(p1, pb1, pb2, p2);
}

namespace
{
    const float kCircleOffset = 4.0f / 3.0f * (sqrtf(2.0f) - 1.0f);

    const Vec4f kQuarterB0(1.0f, 1.0f, kCircleOffset, 0.0f);
    const Vec4f kQuarterB1(0.0f, kCircleOffset, 1.0f, 1.0f);
}

cSpline2 SplineLib::QuadrantSpline(const Vec2f& p, float r, int quadrant)
{
    SL_ASSERT(quadrant >= 0 && quadrant < 4);
    cSpline2 s;

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

    s.xb = r * s.xb + Vec4f(p.x, p.x, p.x, p.x);
    s.yb = r * s.yb + Vec4f(p.y, p.y, p.y, p.y);

    return s;
}

void SplineLib::CircleSplines(const Vec2f& p, float r, cSpline2 splines[4])
{
    for (int i = 0; i < 4; i++)
        splines[i] = QuadrantSpline(p, r, i);
}

namespace
{
    inline cSpline2 SplineFromPoints2(const char* p8, size_t stride, int i0, int i1, int i2, int i3, float tension)
    {
        Vec2f p0 = *(Vec2f*) (p8 + i0 * stride);
        Vec2f p1 = *(Vec2f*) (p8 + i1 * stride);
        Vec2f p2 = *(Vec2f*) (p8 + i2 * stride);
        Vec2f p3 = *(Vec2f*) (p8 + i3 * stride);

        float s = (1.0f - tension) * (1.0f / 6.0f);

        Vec2f pb1 = p1 + s * (p2 - p0);
        Vec2f pb2 = p2 - s * (p3 - p1);

        return BezierSpline(p1, pb1, pb2, p2);
    }
}

int SplineLib::SplinesFromPoints(int numPoints, const Vec2f pi[], int maxSplines, cSpline2 splines[], float tension, size_t stride)
{
    SL_ASSERT(numPoints >= 0);
    SL_ASSERT(maxSplines >= 0 && maxSplines >= NumSplinesForPoints(numPoints));

    const char* p8 = (const char*) pi;

    switch (numPoints)
    {
    case 0:
        return 0;
    case 1:
        *splines = SplineFromPoints2(p8, stride, 0, 0, 0, 0, tension);
        return 1;
    case 2:
        *splines = SplineFromPoints2(p8, stride, 0, 0, 1, 1, tension);
        return 1;
    }

    *splines++ = SplineFromPoints2(p8, stride, 0, 0, 1, 2, tension);

    for (int i = 0; i < numPoints - 3; i++)
    {
        *splines++ = SplineFromPoints2(p8, stride, 0, 1, 2, 3, tension);
        p8 += stride;
    }

    *splines++ = SplineFromPoints2(p8, stride, 0, 1, 2, 2, tension);

    return numPoints - 1;
}

int SplineLib::SplinesFromBezier(int numPoints, const Vec2f points[], const Vec2f hullPoints[], cSpline2 splines[], bool split)
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

int SplineLib::SplinesFromHermite(int numPoints, const Vec2f points[], const Vec2f tangents  [], cSpline2 splines[], bool split)
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

namespace
{
    inline Vec2f Evaluate(const cSpline2& spline, const Vec4f& w)
    /// Evaluate spline with given weights
    {
        return Vec2f
        (
            dot(spline.xb, w),
            dot(spline.yb, w)
        );
    }
}

Vec2f SplineLib::Position(const cSpline2& spline, float t)
{
    return Evaluate(spline, BezierWeights(t));
}

Vec2f SplineLib::Velocity(const cSpline2& spline, float t)
{
    Vec4f dt4(0, 1, 2 * t, 3 * t * t);

    return Evaluate(spline, BezierWeights(dt4));
}

Vec2f SplineLib::Acceleration(const cSpline2& spline, float t)
{
    Vec4f ddt4(0, 0, 2, 6 * t);

    return Evaluate(spline, BezierWeights(ddt4));
}

float SplineLib::Curvature(const cSpline2& spline, float t)
{
    Vec2f v = Velocity    (spline, t);
    Vec2f a = Acceleration(spline, t);

    float avCrossLen = fabsf(v.x * a.y - v.y * a.x);
    float vLen = len(v);

    if (vLen == 0.0f)
        return 1e10f;

    return avCrossLen / (vLen * vLen * vLen);
}

void SplineLib::Frame(const cSpline2& spline, float t, Mat2f* frameOut)
{
    Vec2f v = Velocity    (spline, t);

    Mat2f& frame = *frameOut;

    frame.rows[0] = norm_safe(v);
    frame.rows[1] = cross(frame.rows[0]);
}

float SplineLib::LengthEstimate(const cSpline2& s, float* error)
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

float SplineLib::Length(const cSpline2& s, float maxError)
{
    float error;
    float length = LengthEstimate(s, &error);

    if (error > maxError)
    {
        cSpline2 s0;
        cSpline2 s1;

        Split(s, &s0, &s1);

        return Length(s0, maxError) + Length(s1, maxError);
    }

    return length;
}

float SplineLib::Length(const cSpline2& s, float t0, float t1, float maxError)
{
    SL_ASSERT(t0 >= 0.0f && t0 <  1.0f);
    SL_ASSERT(t1 >= 0.0f && t1 <= 1.0f);
    SL_ASSERT(t0 <= t1);

    cSpline2 s0, s1;

    if (t0 == 0.0f)
    {
        if (t1 == 1.0f)
            return Length(s, maxError);

        Split(s, t1, &s0, &s1);
        return Length(s0, maxError);
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

Bounds2f SplineLib::FastBounds(const cSpline2& spline)
{
    Vec2f bx = HullBounds(spline.xb);
    Vec2f by = HullBounds(spline.yb);

    Bounds2f result = { { bx.x, by.x }, { bx.y, by.y } };
    return result;
}

Bounds2f SplineLib::ExactBounds(const cSpline2& spline)
{
    Vec2f bx = ::ExactBounds(spline.xb);
    Vec2f by = ::ExactBounds(spline.yb);

    Bounds2f result = { { bx.x, by.x }, { bx.y, by.y } };
    return result;
}


void SplineLib::Split(const cSpline2& spline, cSpline2* spline0, cSpline2* spline1)
{
    ::Split(spline.xb, &spline0->xb, &spline1->xb);
    ::Split(spline.yb, &spline0->yb, &spline1->yb);
}

void SplineLib::Split(const cSpline2& spline, float t, cSpline2* spline0, cSpline2* spline1)
{
    ::Split(t, spline.xb, &spline0->xb, &spline1->xb);
    ::Split(t, spline.yb, &spline0->yb, &spline1->yb);
}

bool SplineLib::Join(const cSpline2& s0, const cSpline2& s1, cSpline2* splineOut)
{
    return
       ::Join(s0.xb, s1.xb, &splineOut->xb)
    && ::Join(s0.yb, s1.yb, &splineOut->yb);
}

void SplineLib::Split(vector<cSpline2>* splinesIn)
{
    vector<cSpline2> splines;

    for (const cSpline2& s : *splinesIn)
    {
        cSpline2 s0, s1;

        Split(s, &s0, &s1);
        splines.push_back(s0);
        splines.push_back(s1);
    }

    splinesIn->swap(splines);
}

void SplineLib::Split(vector<cSpline2>* splinesIn, int n)
{
    vector<cSpline2> splines;

    for (const cSpline2& s : *splinesIn)
    {
        cSpline2 ss(s);
        cSpline2 s0, s1;

        for (int i = n; i > 1; i--)
        {
            Split(ss, 1.0f / i, &s0, &ss);
            splines.push_back(s0);
        }
        splines.push_back(ss);
    }

    splinesIn->swap(splines);
}

void SplineLib::Join(vector<cSpline2>* splinesIn)
{
    vector<cSpline2> splines;
    const cSpline2* prevS = 0;

    for (const cSpline2& s : *splinesIn)
    {
        if (!prevS)
        {
            prevS = &s;
            continue;
        }

        cSpline2 sj;
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

namespace
{
    void SubdivideForLength(const cSpline2& s, vector<cSpline2>* splines, float tolerance)
    {
        float error;
        float length = LengthEstimate(s, &error);

        if (error <= tolerance * length)
            splines->push_back(s);
        else
        {
            cSpline2 s1, s2;
            Split(s, &s1, &s2);

            SubdivideForLength(s1, splines, tolerance);
            SubdivideForLength(s2, splines, tolerance);
        }
    }
}

void SplineLib::SubdivideForLength(vector<cSpline2>* splinesIn, float tolerance)
{
    vector<cSpline2> splines;

    for (const cSpline2& s : *splinesIn)
        ::SubdivideForLength(s, &splines, tolerance);

    splinesIn->swap(splines);
}

namespace
{
    inline float ArcError(const cSpline2& s, float* tSplit)
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

    void SubdivideForT(const cSpline2& s, vector<cSpline2>* splines, float tolerance)
    {
        float splitT;
        float err = ArcError(s, &splitT);

        if (err <= tolerance)
            splines->push_back(s);
        else
        {
            cSpline2 s1, s2;
            Split(s, splitT, &s1, &s2);

            SubdivideForT(s1, splines, tolerance);
            SubdivideForT(s2, splines, tolerance);
        }
    }
}

void SplineLib::SubdivideForT(vector<cSpline2>* splinesIn, float tolerance)
{
    vector<cSpline2> splines;

    for (const cSpline2& s : *splinesIn)
        ::SubdivideForT(s, &splines, tolerance);

    splinesIn->swap(splines);
}

namespace
{
    inline float ClosestPoint(const Vec2f& p, const Vec2f& p0, const Vec2f& p1)
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

    void FindClosestPointNewtonRaphson(const cSpline2& spline, Vec2f p, float sIn, int maxIterations, float* tOut, float* dOut)
    {
        SL_ASSERT(sIn >= 0.0f && sIn <= 1.0f);
        const float maxS = 1.0f - 1e-6f;

        float skLast = sIn;
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

            skLast = sk;
        }

        (*tOut) = sk;
        (*dOut) = dk;
    }
}

float SplineLib::FindClosestPoint(const Vec2f& p, const cSpline2& spline)
{
    // Approximate s from straight line between the start and end.
    float s = ClosestPoint(p, Position0(spline), Position1(spline));

    // Use this as starting point for Newton-Raphson solve.
    float d;
    FindClosestPointNewtonRaphson(spline, p, s, 8, &s, &d);

    return s;
}

float SplineLib::FindClosestPoint(const Vec2f& p, int numSplines, const cSpline2 splines[], int* index)
{
    vector<cSubSpline2> nearbyInfo;

    FindNearbySplines(p, numSplines, splines, &nearbyInfo);
    return FindClosestPoint(p, numSplines, splines, nearbyInfo, index);
}


namespace
{
    void FindMinMaxDistance2s(const Vec2f& p, const Bounds2f& bbox, float* minD2, float* maxD2)
    {
        const Vec2f& p0 = bbox.min;
        const Vec2f& p1 = bbox.max;

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

    void FindMinMaxDistance2s(const Vec2f& p, const cSpline2& spline, float* minD2, float* maxD2)
    {
        Bounds2f bbox = FastBounds(spline);

        FindMinMaxDistance2s(p, bbox, minD2, maxD2);
    }

    void Split(const cSubSpline2& s, cSubSpline2* s0, cSubSpline2* s1)
    {
        ::Split(s.mSpline.xb, &s0->mSpline.xb, &s1->mSpline.xb);
        ::Split(s.mSpline.yb, &s0->mSpline.yb, &s1->mSpline.yb);

        s0->mParent = s.mParent;
        s1->mParent = s.mParent;
    }
}

int SplineLib::FindNearbySplines(const Vec2f& p, int numSplines, const cSpline2 splines[], vector<cSubSpline2>* results, float* smallestFarOut, int numIter)
{
    vector<cSubSpline2>& nearSplines = *results;

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
                    if (nearSplines[j].mD2 < smallestFar)
                    {
                        if (dj < j)
                            nearSplines[dj] = nearSplines[j];

                        dj++;
                    }

                nearSplines.resize(dj);
            }

            cSubSpline2 ss = { splines[i], i, near };
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
                    nearSplines  [di] = nearSplines  [i];

                di++;
            }

        nearSplines.resize(di);
        numNearSplines = di;
    }

    if (smallestFarOut)
        *smallestFarOut = smallestFar;

    return numNearSplines;
}

float SplineLib::FindClosestPoint(const Vec2f& p, int numSplines, const cSpline2 splines[], const vector<cSubSpline2>& nearbySplines, int* index)
{
    int prevParent = -1;
    float minD = FLT_MAX;
    float minT = 0.0f;

    *index = -1;

    for (const cSubSpline2& subSpline : nearbySplines)
    {
        if (subSpline.mParent != prevParent)
        {
            SL_ASSERT(subSpline.mParent >= 0 && subSpline.mParent < numSplines);

            const cSpline2& spline = splines[subSpline.mParent];

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
    struct cSubSplineT2
    {
        cSpline2 mSpline;
        float    mT0;
        float    mT1;
    };

    inline void Split(const cSubSplineT2& s, cSubSplineT2* s0, cSubSplineT2* s1)
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
        const cSubSplineT2& spline0,
        const cSubSplineT2& spline1,
        int   dest,
        int   maxDest,
        float results[][2],
        float tolerance
    )
    {
        SL_ASSERT(dest < maxDest);

        Bounds2f bbox0 = ExactBounds(spline0.mSpline);
        Bounds2f bbox1 = ExactBounds(spline1.mSpline);

        if (!Intersects(bbox0, bbox1))
            return dest;

        if (Larger(bbox0, tolerance))
        {
            cSubSplineT2 spline00, spline01;
            Split(spline0, &spline00, &spline01);

            if (Larger(bbox1, tolerance))
            {
                cSubSplineT2 spline10, spline11;
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
            cSubSplineT2 spline10, spline11;
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

int SplineLib::FindSplineIntersections(const cSpline2& spline0, const cSpline2& spline1, int maxResults, float results[][2], float tolerance)
{
    if (maxResults <= 0)
        return 0;

    cSubSplineT2 subSpline0 = { spline0, 0.0f, 1.0f };
    cSubSplineT2 subSpline1 = { spline1, 0.0f, 1.0f };

    return FindSubSplineIntersections(subSpline0, subSpline1, 0, maxResults, results, tolerance);
}

namespace
{
    int FindSplineIntersections
    (
        int   is0, int numSplines0, const cSpline2 splines0[],
        int   is1, int numSplines1, const cSpline2 splines1[],
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
            Bounds2f b0 = FastBounds(splines0[0]);
            Bounds2f b1 = FastBounds(splines1[0]);

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

        cSubSplineT2 st0, st1;

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

int SplineLib::FindSplineIntersections
(
    int   numSplines0, const cSpline2 splines0[],
    int   numSplines1, const cSpline2 splines1[],
    int   maxResults,
    int   resultsI[][2],
    float resultsT[][2],
    float tolerance
)
{
    return ::FindSplineIntersections(0, numSplines0, splines0, 0, numSplines1, splines1, maxResults, resultsI, resultsT, tolerance);
}

int SplineLib::FindSplineIntersections
(
    int   numSplines, const cSpline2 splines[],
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

    cSubSplineT2 st0, st1;

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

float SplineLib::AdvanceAgent(const cSpline2& spline, float t, float ds)
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

bool SplineLib::AdvanceAgent(int* index, float* t, int numSplines, const cSpline2 splines[], float ds)
{
    *t = AdvanceAgent(splines[*index], *t, ds);

    return ::AdvanceAgent(index, t, numSplines);
}

cSpline2 SplineLib::Reverse(const cSpline2& spline)
{
    return
    {
        Vec4f(spline.xb.w, spline.xb.z, spline.xb.y, spline.xb.x),
        Vec4f(spline.yb.w, spline.yb.z, spline.yb.y, spline.yb.x)
    };
}

void SplineLib::Reverse(vector<cSpline2>* splines)
{
    int n = int(splines->size());
    int h = n / 2;

    for (int i = 0; i < h; i++)
    {
        cSpline2& s0 = (*splines)[i];
        cSpline2& s1 = (*splines)[n - i - 1];

        cSpline2 sr0 = Reverse(s1);
        cSpline2 sr1 = Reverse(s0);

        s0 = sr0;
        s1 = sr1;
    }

    if (2 * h < n)
        (*splines)[h] = Reverse((*splines)[h]);
}

cSpline2 SplineLib::Offset(const cSpline2& spline, float offset)
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

void SplineLib::Offset(vector<cSpline2>* splines, float offset)
{
    for (cSpline2& s : *splines)
        s = Offset(s, offset);
}








////////////////////////////////////////////////////////////////////////////////
// 3D
////////////////////////////////////////////////////////////////////////////////


cSpline3 SplineLib::BezierSpline(const Vec3f& p0, const Vec3f& p1, const Vec3f& p2, const Vec3f& p3)
{
    return cSpline3
    {
        Vec4f(p0.x, p1.x, p2.x, p3.x),
        Vec4f(p0.y, p1.y, p2.y, p3.y),
        Vec4f(p0.z, p1.z, p2.z, p3.z)
    };
}

cSpline3 SplineLib::HermiteSpline(const Vec3f& p0, const Vec3f& p1, const Vec3f& v0, const Vec3f& v1)
{
    Vec3f pb1 = p0 + (1.0f / 3.0f) * v0;
    Vec3f pb2 = p1 - (1.0f / 3.0f) * v1;

    return BezierSpline(p0, pb1, pb2, p1);
}

cSpline3 SplineLib::CatmullRomSpline(const Vec3f& p0, const Vec3f& p1, const Vec3f& p2, const Vec3f& p3)
{
    Vec3f pb1 = p1 + (1.0f / 6.0f) * (p2 - p0);
    Vec3f pb2 = p2 - (1.0f / 6.0f) * (p3 - p1);

    return BezierSpline(p1, pb1, pb2, p2);
}

namespace
{
    inline cSpline3 SplineFromPoints3(const char* p8, size_t stride, int i0, int i1, int i2, int i3, float tension)
    {
        Vec3f p0 = *(Vec3f*) (p8 + i0 * stride);
        Vec3f p1 = *(Vec3f*) (p8 + i1 * stride);
        Vec3f p2 = *(Vec3f*) (p8 + i2 * stride);
        Vec3f p3 = *(Vec3f*) (p8 + i3 * stride);

        float s = (1.0f - tension) * (1.0f / 6.0f);

        Vec3f pb1 = p1 + s * (p2 - p0);
        Vec3f pb2 = p2 - s * (p3 - p1);

        return BezierSpline(p1, pb1, pb2, p2);
    }
}

int SplineLib::SplinesFromPoints(int numPoints, const Vec3f pi[], int maxSplines, cSpline3 splines[], float tension, size_t stride)
{
    SL_ASSERT(numPoints >= 0);
    SL_ASSERT(maxSplines >= 0 && maxSplines >= NumSplinesForPoints(numPoints));

    const char* p8 = (const char*) pi;

    switch (numPoints)
    {
    case 0:
        return 0;
    case 1:
        *splines = SplineFromPoints3(p8, stride, 0, 0, 0, 0, tension);
        return 1;
    case 2:
        *splines = SplineFromPoints3(p8, stride, 0, 0, 1, 1, tension);
        return 1;
    }

    *splines++ = SplineFromPoints3(p8, stride, 0, 0, 1, 2, tension);

    for (int i = 0; i < numPoints - 3; i++)
    {
        *splines++ = SplineFromPoints3(p8, stride, 0, 1, 2, 3, tension);
        p8 += stride;
    }

    *splines++ = SplineFromPoints3(p8, stride, 0, 1, 2, 2, tension);

    return numPoints - 1;
}

int SplineLib::SplinesFromBezier(int numPoints, const Vec3f points[], const Vec3f hullPoints[], cSpline3 splines[], bool split)
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

int SplineLib::SplinesFromHermite(int numPoints, const Vec3f points[], const Vec3f tangents  [], cSpline3 splines[], bool split)
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

namespace
{
    inline Vec3f Evaluate(const cSpline3& spline, const Vec4f& w)
    /// Evaluate spline with given weights
    {
        return Vec3f
        (
            dot(spline.xb, w),
            dot(spline.yb, w),
            dot(spline.zb, w)
        );
    }
}

Vec3f SplineLib::Position(const cSpline3& spline, float t)
{
    return Evaluate(spline, BezierWeights(t));
}

Vec3f SplineLib::Velocity(const cSpline3& spline, float t)
{
    Vec4f dt4(0, 1, 2 * t, 3 * t * t);

    return Evaluate(spline, BezierWeights(dt4));
}

Vec3f SplineLib::Acceleration(const cSpline3& spline, float t)
{
    Vec4f ddt4(0, 0, 2, 6 * t);

    return Evaluate(spline, BezierWeights(ddt4));
}

float SplineLib::Curvature(const cSpline3& spline, float t)
{
    Vec3f v = Velocity    (spline, t);
    Vec3f a = Acceleration(spline, t);

    float avCrossLen = len(cross(v, a));
    float vLen = len(v);

    if (vLen == 0.0f)
        return 1e10f;

    return avCrossLen / (vLen * vLen * vLen);
}

void SplineLib::Frame(const cSpline3& spline, float t, Mat3f* frameOut)
{
    Vec3f v = Velocity    (spline, t);
    Vec3f a = Acceleration(spline, t);

    Vec3f g = norm_safe(v);             // tangent
    Vec3f b = norm_safe(cross(v, a));
    Vec3f n = cross(g, b);

    Mat3f& frame = *frameOut;

    frame.rows[0] = n;   // x = right
    frame.rows[1] = g;   // y = forward
    frame.rows[2] = b;   // z = up = normal
}

void SplineLib::FrameSmooth(const cSpline3& spline, float t, Mat3f* frameOut)
{
    // Issues:
    //   'b' can flip sign as 'a' moves from one side of the curve to the other
    //   'a' can go to zero on flat sections of the curve. (In which case we should maintain previous values.)
    //   'v' can go to zero at inflection points (ditto).
    // Generally we want to preserve the tangent direction above everything, then the up vector, and last the right vector.
    // E.g., if agent is going into a point and coming back out, we desire that the right vector flips.

    Vec3f v = Velocity    (spline, t);
    Vec3f a = Acceleration(spline, t);

    Mat3f& frame = *frameOut;

    Vec3f g;
    if (sqrlen(v) >= 1e-3f)
        g = norm(v);
    else
        g = frame.rows[0];   // fall back to keeping existing tangent

    Vec3f ga = cross(g, a);
    Vec3f b;
    if (sqrlen(ga) >= 1e-3f)
    {
        b = norm(ga);

        if (dot(b, frame.rows[2]) < 0.0f)
            b = -b; // correct b
    }
    else
        b = frame.rows[2];   // fall back to previous up

    Vec3f n = cross(g, b);

    frame.rows[0] = n;   // x = right
    frame.rows[1] = g;   // y = forward
    frame.rows[2] = b;   // z = up = normal
}


void SplineLib::FrameUp(const cSpline3& spline, float t, const Vec3f& up, Mat3f* frameOut)
{
    Vec3f v = Velocity(spline, t);

    Vec3f g = norm_safe(v);             // tangent
    Vec3f n = norm_safe(cross(g, up));
    Vec3f b = cross(n, g);

    Mat3f& frame = *frameOut;

    frame.rows[0] = n;   // x = right
    frame.rows[1] = g;   // y = forward
    frame.rows[2] = b;   // z = up = normal
}

void SplineLib::FrameZUp(const cSpline3& spline, float t, Mat3f* frameOut)
{
    Vec3f v = Velocity(spline, t);

    Vec3f g = norm_safe(v);             // tangent
    Vec3f n = norm_safe(cross_z(g));
    Vec3f b = cross(n, g);

    Mat3f& frame = *frameOut;

    frame.rows[0] = n;   // x = right
    frame.rows[1] = g;   // y = forward
    frame.rows[2] = b;   // z = up = normal
}

float SplineLib::LengthEstimate(const cSpline3& s, float* error)
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

float SplineLib::Length(const cSpline3& s, float maxError)
{
    float error;
    float length = LengthEstimate(s, &error);

    if (error > maxError)
    {
        cSpline3 s0;
        cSpline3 s1;

        Split(s, &s0, &s1);

        return Length(s0, maxError) + Length(s1, maxError);
    }

    return length;
}

float SplineLib::Length(const cSpline3& s, float t0, float t1, float maxError)
{
    SL_ASSERT(t0 >= 0.0f && t0 <  1.0f);
    SL_ASSERT(t1 >= 0.0f && t1 <= 1.0f);
    SL_ASSERT(t0 <= t1);

    cSpline3 s0, s1;

    if (t0 == 0.0f)
    {
        if (t1 == 1.0f)
            return Length(s, maxError);

        Split(s, t1, &s0, &s1);
        return Length(s0, maxError);
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

Bounds3f SplineLib::FastBounds(const cSpline3& spline)
{
    Vec2f bx = HullBounds(spline.xb);
    Vec2f by = HullBounds(spline.yb);
    Vec2f bz = HullBounds(spline.zb);

    Bounds3f result = { { bx.x, by.x, bz.x }, { bx.y, by.y, bz.y } };
    return result;
}

Bounds3f SplineLib::ExactBounds(const cSpline3& spline)
{
    Vec2f bx = ::ExactBounds(spline.xb);
    Vec2f by = ::ExactBounds(spline.yb);
    Vec2f bz = ::ExactBounds(spline.zb);

    Bounds3f result = { { bx.x, by.x, bz.x }, { bx.y, by.y, bz.y } };
    return result;
}


void SplineLib::Split(const cSpline3& spline, cSpline3* spline0, cSpline3* spline1)
{
    ::Split(spline.xb, &spline0->xb, &spline1->xb);
    ::Split(spline.yb, &spline0->yb, &spline1->yb);
    ::Split(spline.zb, &spline0->zb, &spline1->zb);
}

void SplineLib::Split(const cSpline3& spline, float t, cSpline3* spline0, cSpline3* spline1)
{
    ::Split(t, spline.xb, &spline0->xb, &spline1->xb);
    ::Split(t, spline.yb, &spline0->yb, &spline1->yb);
    ::Split(t, spline.zb, &spline0->zb, &spline1->zb);
}

bool SplineLib::Join(const cSpline3& s0, const cSpline3& s1, cSpline3* splineOut)
{
    return
       ::Join(s0.xb, s1.xb, &splineOut->xb)
    && ::Join(s0.yb, s1.yb, &splineOut->yb)
    && ::Join(s0.zb, s1.zb, &splineOut->zb);
}

void SplineLib::Split(vector<cSpline3>* splinesIn)
{
    vector<cSpline3> splines;

    for (const cSpline3& s : *splinesIn)
    {
        cSpline3 s0, s1;

        Split(s, &s0, &s1);
        splines.push_back(s0);
        splines.push_back(s1);
    }

    splinesIn->swap(splines);
}

void SplineLib::Split(vector<cSpline3>* splinesIn, int n)
{
    vector<cSpline3> splines;

    for (const cSpline3& s : *splinesIn)
    {
        cSpline3 ss(s);
        cSpline3 s0, s1;

        for (int i = n; i > 1; i--)
        {
            Split(ss, 1.0f / i, &s0, &ss);
            splines.push_back(s0);
        }
        splines.push_back(ss);
    }

    splinesIn->swap(splines);
}

void SplineLib::Join(vector<cSpline3>* splinesIn)
{
    vector<cSpline3> splines;
    const cSpline3* prevS = 0;

    for (const cSpline3& s : *splinesIn)
    {
        if (!prevS)
        {
            prevS = &s;
            continue;
        }

        cSpline3 sj;
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

namespace
{
    void SubdivideForLength(const cSpline3& s, vector<cSpline3>* splines, float tolerance)
    {
        float error;
        float length = LengthEstimate(s, &error);

        if (error <= tolerance * length)
            splines->push_back(s);
        else
        {
            cSpline3 s1, s2;
            Split(s, &s1, &s2);

            SubdivideForLength(s1, splines, tolerance);
            SubdivideForLength(s2, splines, tolerance);
        }
    }
}

void SplineLib::SubdivideForLength(vector<cSpline3>* splinesIn, float tolerance)
{
    vector<cSpline3> splines;

    for (const cSpline3& s : *splinesIn)
        ::SubdivideForLength(s, &splines, tolerance);

    splinesIn->swap(splines);
}

namespace
{
    float ArcError(const cSpline3& s, float* tSplit)
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

    void SubdivideForT(const cSpline3& s, vector<cSpline3>* splines, float tolerance)
    {
        float splitT;
        float err = ArcError(s, &splitT);

        if (err <= tolerance)
            splines->push_back(s);
        else
        {
            cSpline3 s1, s2;
            Split(s, splitT, &s1, &s2);

            SubdivideForT(s1, splines, tolerance);
            SubdivideForT(s2, splines, tolerance);
        }
    }
}

void SplineLib::SubdivideForT(vector<cSpline3>* splinesIn, float tolerance)
{
    vector<cSpline3> splines;

    for (const cSpline3& s : *splinesIn)
        ::SubdivideForT(s, &splines, tolerance);

    splinesIn->swap(splines);
}

namespace
{
    inline float ClosestPoint(const Vec3f& p, const Vec3f& p0, const Vec3f& p1)
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

    void FindClosestPointNewtonRaphson(const cSpline3& spline, Vec3f p, float sIn, int maxIterations, float* tOut, float* dOut)
    {
        SL_ASSERT(sIn >= 0.0f && sIn <= 1.0f);
        const float maxS = 1.0f - 1e-6f;

        float skLast = sIn;
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

            skLast = sk;
        }

        (*tOut) = sk;
        (*dOut) = dk;
    }
}

float SplineLib::FindClosestPoint(const Vec3f& p, const cSpline3& spline)
{
    // Approximate s from straight line between the start and end.
    float s = ClosestPoint(p, Position0(spline), Position1(spline));

    // Use this as starting point for Newton-Raphson solve.
    float d;
    FindClosestPointNewtonRaphson(spline, p, s, 8, &s, &d);

    return s;
}

float SplineLib::FindClosestPoint(const Vec3f& p, int numSplines, const cSpline3 splines[], int* index)
{
    vector<cSubSpline3> nearbyInfo;

    FindNearbySplines(p, numSplines, splines, &nearbyInfo);
    return FindClosestPoint(p, numSplines, splines, nearbyInfo, index);
}


namespace
{
    void FindMinMaxDistance2s(const Vec3f& p, const Bounds3f& bbox, float* minD2, float* maxD2)
    {
        const Vec3f& p0 = bbox.min;
        const Vec3f& p1 = bbox.max;

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

    void FindMinMaxDistance2s(const Vec3f& p, const cSpline3& spline, float* minD2, float* maxD2)
    {
        Bounds3f bbox = FastBounds(spline);

        FindMinMaxDistance2s(p, bbox, minD2, maxD2);
    }

    void Split(const cSubSpline3& s, cSubSpline3* s0, cSubSpline3* s1)
    {
        ::Split(s.mSpline.xb, &s0->mSpline.xb, &s1->mSpline.xb);
        ::Split(s.mSpline.yb, &s0->mSpline.yb, &s1->mSpline.yb);
        ::Split(s.mSpline.zb, &s0->mSpline.zb, &s1->mSpline.zb);

        s0->mParent = s.mParent;
        s1->mParent = s.mParent;
    }
}

int SplineLib::FindNearbySplines(const Vec3f& p, int numSplines, const cSpline3 splines[], vector<cSubSpline3>* results, float* smallestFarOut, int numIter)
{
    vector<cSubSpline3>& nearSplines = *results;

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
                    if (nearSplines[j].mD2 < smallestFar)
                    {
                        if (dj < j)
                            nearSplines[dj] = nearSplines[j];

                        dj++;
                    }

                nearSplines.resize(dj);
            }

            cSubSpline3 ss = { splines[i], i, near };
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
                    nearSplines  [di] = nearSplines  [i];

                di++;
            }

        nearSplines.resize(di);
        numNearSplines = di;
    }

    if (smallestFarOut)
        *smallestFarOut = smallestFar;

    return numNearSplines;
}

float SplineLib::FindClosestPoint(const Vec3f& p, int numSplines, const cSpline3 splines[], const vector<cSubSpline3>& nearbySplines, int* index)
{
    int prevParent = -1;
    float minD = FLT_MAX;
    float minT = 0.0f;

    *index = -1;

    for (const cSubSpline3& subSpline : nearbySplines)
    {
        if (subSpline.mParent != prevParent)
        {
            SL_ASSERT(subSpline.mParent >= 0 && subSpline.mParent < numSplines);

            const cSpline3& spline = splines[subSpline.mParent];

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
    struct cSubSplineT3
    {
        cSpline3 mSpline;
        float    mT0;
        float    mT1;
    };

    inline void Split(const cSubSplineT3& s, cSubSplineT3* s0, cSubSplineT3* s1)
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
        const cSubSplineT3& spline0,
        const cSubSplineT3& spline1,
        int   dest,
        int   maxDest,
        float results[][2],
        float tolerance
    )
    {
        SL_ASSERT(dest < maxDest);

        Bounds3f bbox0 = ExactBounds(spline0.mSpline);
        Bounds3f bbox1 = ExactBounds(spline1.mSpline);

        if (!Intersects(bbox0, bbox1))
            return dest;

        if (Larger(bbox0, tolerance))
        {
            cSubSplineT3 spline00, spline01;
            Split(spline0, &spline00, &spline01);

            if (Larger(bbox1, tolerance))
            {
                cSubSplineT3 spline10, spline11;
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
            cSubSplineT3 spline10, spline11;
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

int SplineLib::FindSplineIntersections(const cSpline3& spline0, const cSpline3& spline1, int maxResults, float results[][2], float tolerance)
{
    if (maxResults <= 0)
        return 0;

    cSubSplineT3 subSpline0 = { spline0, 0.0f, 1.0f };
    cSubSplineT3 subSpline1 = { spline1, 0.0f, 1.0f };

    return FindSubSplineIntersections(subSpline0, subSpline1, 0, maxResults, results, tolerance);
}

namespace
{
    int FindSplineIntersections
    (
        int   is0, int numSplines0, const cSpline3 splines0[],
        int   is1, int numSplines1, const cSpline3 splines1[],
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
            Bounds3f b0 = FastBounds(splines0[0]);
            Bounds3f b1 = FastBounds(splines1[0]);

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

        cSubSplineT3 st0, st1;

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

int SplineLib::FindSplineIntersections
(
    int   numSplines0, const cSpline3 splines0[],
    int   numSplines1, const cSpline3 splines1[],
    int   maxResults,
    int   resultsI[][2],
    float resultsT[][2],
    float tolerance
)
{
    return ::FindSplineIntersections(0, numSplines0, splines0, 0, numSplines1, splines1, maxResults, resultsI, resultsT, tolerance);
}

int SplineLib::FindSplineIntersections
(
    int   numSplines, const cSpline3 splines[],
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

    cSubSplineT3 st0, st1;

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

float SplineLib::AdvanceAgent(const cSpline3& spline, float t, float ds)
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

bool SplineLib::AdvanceAgent(int* index, float* t, int numSplines, const cSpline3 splines[], float ds)
{
    *t = AdvanceAgent(splines[*index], *t, ds);

    return ::AdvanceAgent(index, t, numSplines);
}

cSpline3 SplineLib::Reverse(const cSpline3& spline)
{
    return
    {
        Vec4f(spline.xb.w, spline.xb.z, spline.xb.y, spline.xb.x),
        Vec4f(spline.yb.w, spline.yb.z, spline.yb.y, spline.yb.x),
        Vec4f(spline.zb.w, spline.zb.z, spline.zb.y, spline.zb.x)
    };
}

void SplineLib::Reverse(vector<cSpline3>* splines)
{
    int n = int(splines->size());
    int h = n / 2;

    for (int i = 0; i < h; i++)
    {
        cSpline3& s0 = (*splines)[i];
        cSpline3& s1 = (*splines)[n - i - 1];

        cSpline3 sr0 = Reverse(s1);
        cSpline3 sr1 = Reverse(s0);

        s0 = sr0;
        s1 = sr1;
    }

    if (2 * h < n)
        (*splines)[h] = Reverse((*splines)[h]);
}

cSpline3 SplineLib::Offset(const cSpline3& spline, float offset)
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

void SplineLib::Offset(vector<cSpline3>* splines, float offset)
{
    for (cSpline3& s : *splines)
        s = Offset(s, offset);
}
