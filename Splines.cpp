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
    // Mini VL lib
    inline float sqr(float x) { return x * x; }
    inline float lerp(float a, float b, float t) { return (1.0f - t) * a + t * b; }

    inline float dot(Vec2f a, Vec2f b)       { return a.x * b.x + a.y * b.y; }
    inline float len   (Vec2f v)             { return sqrtf(v.x * v.x + v.y * v.y); }
    inline float sqrlen(Vec2f v)             { return v.x * v.x + v.y * v.y; }
    inline Vec2f abs   (Vec2f v)             { return { fabsf(v.x), fabsf(v.y) }; }

    inline Vec2f operator+(Vec2f a, Vec2f b) { return { a.x + b.x, a.y + b.y }; }
    inline Vec2f operator-(Vec2f a, Vec2f b) { return { a.x - b.x, a.y - b.y }; }
    inline Vec2f operator*(float s, Vec2f a) { return { s * a.x, s * a.y }; }

    inline float dot   (Vec3f a, Vec3f b)    { return a.x * b.x + a.y * b.y + a.z * b.z; }
    inline float len   (Vec3f v)             { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); }
    inline float sqrlen(Vec3f v)             { return v.x * v.x + v.y * v.y + v.z * v.z; }
    inline Vec3f abs   (Vec3f v)             { return { fabsf(v.x), fabsf(v.y), fabsf(v.z) }; }
    inline Vec3f cross (Vec3f a, Vec3f b)    { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }

    inline Vec3f operator+(Vec3f a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z}; }
    inline Vec3f operator-(Vec3f a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z}; }
    inline Vec3f operator*(float s, Vec3f a) { return { s   * a.x, s   * a.y, s   * a.z}; }

    inline float dot   (Vec4f a, Vec4f b)    { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; } 

    inline float Max(float a, float b)
    {
        return b < a ? a : b;
    }

    inline float Min(float a, float b)
    {
        return a < b ? a : b;
    }

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

    return avCrossLen / (vLen * vLen * vLen);
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
        
        // TODO: strictly speaking should be able to use knot info to better estimate where to split?
        Split(s, &s0, &s1);
        
        return Length(s0, maxError) + Length(s1, maxError);
    }
    
    return length;
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
    
    // Use Newton-Raphson to root-find
    void FindClosestPointNewtonRaphson(const cSpline2& spline, Vec2f p, float sIn, int maxIterations, float* tOut, float* dOut)
    {
        SL_ASSERT(sIn >= 0.0f && sIn <= 1.0f);
        
        float maxS = 1.0f - 1e-6f;
        
        // Newton iteration.
        float skLast = sIn;
        float sk = sIn;
        
        float dk = len(Position(spline, sk) - p);
        
        constexpr float width    = 1e-3f;

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
    std::vector<cSubSpline2> nearbyInfo;
    
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
}

int SplineLib::FindNearbySplines(const Vec2f& p, int numSplines, const cSpline2 splines[], std::vector<cSubSpline2>* results, float* smallestFarOut, int numIter)
{
    std::vector<cSubSpline2>& nearSplines = *results;
    
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

            cSubSpline2 ss = { splines[i], near, i };
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
        {
            Split(nearSplines[i].mSpline, &nearSplines[2 * i].mSpline, &nearSplines[2 * i + 1].mSpline);
            
            nearSplines[2 * i + 0].mParent = nearSplines[i].mParent; 
            nearSplines[2 * i + 1].mParent = nearSplines[i].mParent; 
        }

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

float SplineLib::FindClosestPoint(const Vec2f& p, int numSplines, const cSpline2 splines[], const std::vector<cSubSpline2>& nearbySplines, int* index)
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


////////////////////////////////////////////////////////////////////////////////
// 3D

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

    return avCrossLen / (vLen * vLen * vLen);
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
        // Split and recurse!
        cSpline3 s0;
        cSpline3 s1;

        // TODO: strictly speaking should be able to use knot info to better estimate where to split?
        Split(s, &s0, &s1);

        return Length(s0, maxError) + Length(s1, maxError);
    }

    return length;
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
    && ::Join(s0.yb, s1.yb, &splineOut->yb);
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
        float skLast = sIn;
        float sk = sIn;

        float dk = len(Position(spline, sk) - p);

        constexpr float width    = 1e-3f;
        
        float maxJump  = 0.5f;   // avoid jumping too far, leads to oscillation

        for (int i = 0; i < maxIterations; i++)
        {
            float ss = Clamp(sk, width, 1.0f - width); // so can interpolate points for Newtons method

            float d1 = len(Position(spline, ss - width) - p);
            float d2 = len(Position(spline, ss        ) - p);
            float d3 = len(Position(spline, ss + width) - p);

            // TODO: must be able to do this analytically if assume same spline.
            float g1 = (d2 - d1) / width;
            float g2 = (d3 - d2) / width;

            float grad = (d3 - d1) / (2.0f * width);
            float curv = (g2 - g1) / width;

            float sn;

            if (curv > 0.0f)    // if d' is heading towards a minima, apply NR for d'
                sn = ss - grad / curv; 
            else if (grad != 0.0f)
                sn = ss - d2 / grad; // otherwise, apply for D.
            else
                sn = sk;

            sn = Clamp(sn, sk - maxJump, sk + maxJump);   // avoid large steps, often unstable.

            // only update our estimate if the new value is in range and closer.
            if (sn >= 0.0f && sn < 1.0f)
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

float SplineLib::FindClosestPoint(const cSpline3& spline, const Vec3f& p)
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
    std::vector<cSubSpline3> nearbyInfo;
    
    FindNearbySplines(p, numSplines, splines, &nearbyInfo);
    return FindClosestPoint(p, numSplines, splines, nearbyInfo, index);
}

    
namespace
{
    void FindMinMaxDistance2s(const Vec3f& p, const Bounds3f& bbox, float* minD2, float* maxD2)
    {
        // Find the nearest point to p inside the bbox
        // This can be a bbox vertex, a point on an edge or face, or p itself if it's inside the box
        float minX = Clamp(p.x, bbox.min.x, bbox.max.x);
        float minY = Clamp(p.y, bbox.min.y, bbox.max.y);
        float minZ = Clamp(p.z, bbox.min.z, bbox.max.z);

        // Find the farthest point from p inside the bbox
        // This is always a bbox vertex.
        const Vec3f& p0 = bbox.min;
        const Vec3f& p1 = bbox.max;

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
        Bounds3f bbox = ExactBounds(spline);

        FindMinMaxDistance2s(p, bbox, minD2, maxD2);
    }
}

int SplineLib::FindNearbySplines(const Vec3f& p, int numSplines, const cSpline3 splines[], std::vector<cSubSpline3>* results, float* smallestFarOut, int numIter)
{
    std::vector<cSubSpline3>& nearSplines = *results;
    
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

            cSubSpline3 ss = { splines[i], near, i };
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
        {
            Split(nearSplines[i].mSpline, &nearSplines[2 * i].mSpline, &nearSplines[2 * i + 1].mSpline);
            
            nearSplines[2 * i + 0].mParent = nearSplines[i].mParent; 
            nearSplines[2 * i + 1].mParent = nearSplines[i].mParent; 
        }

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

float SplineLib::FindClosestPoint(const Vec3f& p, int numSplines, const cSpline3 splines[], const std::vector<cSubSpline3>& nearbySplines, int* index)
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

