#ifndef VL_MINI_H
#define VL_MINI_H

#include <math.h>

#define VL_ASSERT(X)

template<class T> inline T sqr(T x) { return x * x; }

template<class T> inline T vl_min(T a, T b) { return a < b ? a : b; }
template<class T> inline T vl_max(T a, T b) { return a > b ? a : b; }

constexpr float vlf_pi = 3.14159265358979323846f;
constexpr float vlf_halfPi = 0.5 * vlf_pi;
constexpr float vlf_twoPi = 2 * vlf_pi;

#define VI(T) T operator[](int i) const { return (&x)[i]; } T& operator[](int i){ return (&x)[i]; }

inline float len   (float x) { return fabsf(x); }
inline float sqrlen(float x) { return x * x; }

struct Vec2f { float x; float y;                    Vec2f() {}; Vec2f(float xi, float yi)                     : x(xi), y(yi)               {}; VI(float) };
struct Vec3f { float x; float y; float z;           Vec3f() {}; Vec3f(float xi, float yi, float zi)           : x(xi), y(yi), z(zi)        {}; VI(float) };
struct Vec4f { float x; float y; float z; float w;  Vec4f() {}; Vec4f(float xi, float yi, float zi, float wi) : x(xi), y(yi), z(zi), w(wi) {}; VI(float) };

struct Vec2d { double x; double y;                      Vec2d() {}; Vec2d(double xi, double yi)                       : x(xi), y(yi)               {}; VI(double) };
struct Vec3d { double x; double y; double z;            Vec3d() {}; Vec3d(double xi, double yi, double zi)            : x(xi), y(yi), z(zi)        {}; VI(double) };
struct Vec4d { double x; double y; double z; double w;  Vec4d() {}; Vec4d(double xi, double yi, double zi, double wi) : x(xi), y(yi), z(zi), w(wi) {}; VI(double) };

struct Mat2f { Vec2f x; Vec2f y;                    Mat2f() {}; Mat2f(Vec2f xi, Vec2f yi)                     : x(xi), y(yi)               {}; VI(Vec2f) };
struct Mat3f { Vec3f x; Vec3f y; Vec3f z;           Mat3f() {}; Mat3f(Vec3f xi, Vec3f yi, Vec3f zi)           : x(xi), y(yi), z(zi)        {}; VI(Vec3f) };
struct Mat4f { Vec4f x; Vec4f y; Vec4f z; Vec4f w;  Mat4f() {}; Mat4f(Vec4f xi, Vec4f yi, Vec4f zi, Vec4f wi) : x(xi), y(yi), z(zi), w(wi) {}; VI(Vec4f) };

inline Vec2f  operator+ (Vec2f  v)          { return { +v.x, +v.y }; }
inline Vec2f  operator- (Vec2f  v)          { return { -v.x, -v.y }; }
inline Vec2f  operator+ (Vec2f  a, Vec2f b) { return { a.x + b.x, a.y + b.y }; }
inline Vec2f  operator- (Vec2f  a, Vec2f b) { return { a.x - b.x, a.y - b.y }; }
inline Vec2f  operator* (Vec2f  a, Vec2f b) { return { a.x * b.x, a.y * b.y }; }
inline Vec2f  operator/ (Vec2f  a, Vec2f b) { return { a.x / b.x, a.y / b.y }; }
inline Vec2f  operator* (float  s, Vec2f a) { return { s * a.x, s * a.y }; }
inline Vec2f  operator* (Vec2f  a, float s) { return { s * a.x, s * a.y }; }
inline Vec2f  operator/ (float  s, Vec2f a) { return { s / a.x, s / a.y }; }
inline Vec2f  operator/ (Vec2f  a, float s) { return { s / a.x, s / a.y }; }
inline Vec2f& operator+=(Vec2f& a, Vec2f b) { a.x += b.x; a.y += b.y; return a; }
inline Vec2f& operator-=(Vec2f& a, Vec2f b) { a.x -= b.x; a.y -= b.y; return a; }
inline Vec2f& operator*=(Vec2f& a, Vec2f b) { a.x *= b.x; a.y *= b.y; return a; }
inline Vec2f& operator/=(Vec2f& a, Vec2f b) { a.x /= b.x; a.y /= b.y; return a; }
inline Vec2f& operator*=(Vec2f& a, float s) { a.x *= s  ; a.y *= s  ; return a; }
inline Vec2f& operator/=(Vec2f& a, float s) { a.x /= s  ; a.y /= s  ; return a; }

inline bool   operator==(Vec3f  a, Vec3f b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
inline bool   operator!=(Vec3f  a, Vec3f b) { return a.x != b.x || a.y != b.y || a.z != b.z; }
inline Vec3f  operator+ (Vec3f  v)          { return { +v.x, +v.y, +v.z }; }
inline Vec3f  operator- (Vec3f  v)          { return { -v.x, -v.y, -v.z }; }
inline Vec3f  operator+ (Vec3f  a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline Vec3f  operator- (Vec3f  a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline Vec3f  operator* (Vec3f  a, Vec3f b) { return { a.x * b.x, a.y * b.y, a.z * b.z }; }
inline Vec3f  operator/ (Vec3f  a, Vec3f b) { return { a.x / b.x, a.y / b.y, a.z / b.z }; }
inline Vec3f  operator* (float  s, Vec3f a) { return { s   * a.x, s   * a.y, s   * a.z }; }
inline Vec3f  operator* (Vec3f  a, float s) { return { s   * a.x, s   * a.y, s   * a.z }; }
inline Vec3f  operator/ (float  s, Vec3f a) { return { s   / a.x, s   / a.y, s   / a.z }; }
inline Vec3f  operator/ (Vec3f  a, float s) { return { a.x / s  , a.y / s  , a.z / s   }; }
inline Vec3f& operator+=(Vec3f& a, Vec3f b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
inline Vec3f& operator-=(Vec3f& a, Vec3f b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
inline Vec3f& operator*=(Vec3f& a, Vec3f b) { a.x *= b.x; a.y *= b.y; a.z *= b.z; return a; }
inline Vec3f& operator/=(Vec3f& a, Vec3f b) { a.x /= b.x; a.y /= b.y; a.z /= b.z; return a; }
inline Vec3f& operator*=(Vec3f& a, float s) { a.x *= s  ; a.y *= s  ; a.z *= s  ; return a; }
inline Vec3f& operator/=(Vec3f& a, float s) { a.x /= s  ; a.y /= s  ; a.z /= s  ; return a; }

inline Vec4f  operator+ (Vec4f  v)          { return { +v.x, +v.y, +v.z, +v.w }; }
inline Vec4f  operator- (Vec4f  v)          { return { -v.x, -v.y, -v.z, -v.w }; }
inline Vec4f  operator+ (Vec4f  a, Vec4f b) { return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w }; }
inline Vec4f  operator- (Vec4f  a, Vec4f b) { return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w }; }
inline Vec4f  operator* (Vec4f  a, Vec4f b) { return { a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w }; }
inline Vec4f  operator/ (Vec4f  a, Vec4f b) { return { a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w }; }
inline Vec4f  operator* (float  s, Vec4f a) { return { s   * a.x, s   * a.y, s   * a.z, s   * a.w }; }
inline Vec4f  operator* (Vec4f  a, float s) { return { s   * a.x, s   * a.y, s   * a.z, s   * a.w }; }
inline Vec4f  operator/ (float  s, Vec4f a) { return { s   / a.x, s   / a.y, s   / a.z, s   / a.w }; }
inline Vec4f  operator/ (Vec4f  a, float s) { return { s   / a.x, s   / a.y, s   / a.z, s   / a.w }; }
inline Vec4f& operator+=(Vec4f& a, Vec4f b) { a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w; return a; }
inline Vec4f& operator-=(Vec4f& a, Vec4f b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w; return a; }
inline Vec4f& operator*=(Vec4f& a, Vec4f b) { a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w; return a; }
inline Vec4f& operator/=(Vec4f& a, Vec4f b) { a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w; return a; }
inline Vec4f& operator*=(Vec4f& a, float s) { a.x *= s  ; a.y *= s  ; a.z *= s  ; a.w *= s  ; return a; }
inline Vec4f& operator/=(Vec4f& a, float s) { a.x /= s  ; a.y /= s  ; a.z /= s  ; a.w /= s  ; return a; }

inline float dot      (Vec2f a, Vec2f b) { return a.x * b.x + a.y * b.y; }
inline float len      (Vec2f v)          { return sqrtf(v.x * v.x + v.y * v.y); }
inline float sqrlen   (Vec2f v)          { return       v.x * v.x + v.y * v.y; }
inline Vec2f norm     (Vec2f v)          { return (1.0f / len(v)) * v; }
inline Vec2f norm_safe(Vec2f v)          { return (1.0f / (len(v) + 1e-8f)) * v; }
inline Vec2f abs      (Vec2f v)          { return { fabsf(v.x), fabsf(v.y) }; }
inline Vec2f cross    (Vec2f v)          { return { -v.y, v.x }; }

inline float dot      (Vec3f a, Vec3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline float len      (Vec3f v)          { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); }
inline float sqrlen   (Vec3f v)          { return       v.x * v.x + v.y * v.y + v.z * v.z; }
inline Vec3f norm     (Vec3f v)          { return (1.0f / len(v)) * v; }
inline Vec3f norm_safe(Vec3f v)          { return (1.0f / (len(v) + 1e-8f)) * v; }
inline Vec3f abs      (Vec3f v)          { return { fabsf(v.x), fabsf(v.y), fabsf(v.z) }; }
inline Vec3f cross    (Vec3f a, Vec3f b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
inline Vec3f cross_z  (Vec3f v)          { return { v.y, -v.x, 0.0f }; }

inline float dot      (Vec4f a, Vec4f b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
inline float len      (Vec4f v)          { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w); }
inline float sqrlen   (Vec4f v)          { return       v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w; }
inline Vec4f norm     (Vec4f v)          { return (1.0f / len(v)) * v; }
inline Vec4f norm_safe(Vec4f v)          { return (1.0f / (len(v) + 1e-8f)) * v; }
inline Vec4f abs      (Vec4f v)          { return { fabsf(v.x), fabsf(v.y), fabsf(v.z), fabsf(v.w) }; }

inline float lerp(float a, float b, float t) { return (1.0f - t) * a + t * b; }
inline Vec2f lerp(Vec2f a, Vec2f b, float t) { return (1.0f - t) * a + t * b; }
inline Vec3f lerp(Vec3f a, Vec3f b, float t) { return (1.0f - t) * a + t * b; }
inline Vec4f lerp(Vec4f a, Vec4f b, float t) { return (1.0f - t) * a + t * b; }

// Necessary utilities pulled from https://github.com/andrewwillmott/vl
typedef Vec4f Quatf;

Quatf QuatMult (const Quatf& a, const Quatf& b); // Concatenate quaternions: the result represents applying 'a' then 'b'.
Vec3f QuatApply(const Vec3f& p, const Quatf& q); // Transform point p by applying quaternion q
Quatf QuatConj (const Quatf& q);                 // Quaternion conjugate. if len(q) = 1, this is also the inverse.

Quatf SLerp(const Quatf& q1, const Quatf& q2, float s); // Return spherical interpolation between q1 and q2

Vec3f LogUnit3(const Quatf& q);   // LogUnit variant that omits the last (zero) component
Quatf ExpUnit3(const Vec3f& lq);  // ExpUnit variant that omits the last (zero) component
Vec3f QuatDiff3(const Quatf& a, const Quatf& b);  // Returns LogUnit3(MakeQuat(a, b)) -- useful for rotational velocity, and other log/tangent-space deltas

inline Quatf QuatMult(const Quatf& a, const Quatf& b)
{
    Quatf result;
    result.x = + a.w * b.x + a.z * b.y - a.y * b.z + a.x * b.w;
    result.y = - a.z * b.x + a.w * b.y + a.x * b.z + a.y * b.w;
    result.z = + a.y * b.x - a.x * b.y + a.w * b.z + a.z * b.w;
    result.w = - a.x * b.x - a.y * b.y - a.z * b.z + a.w * b.w;
    return result;
}

inline const Vec3f& xyz(const Quatf& q) { return (const Vec3f&) q; }
inline       Vec3f& xyz(      Quatf& q) { return (      Vec3f&) q; }

inline Vec3f QuatApply(const Vec3f& p, const Quatf& q)
{
    Vec3f b0 = cross(xyz(q), p);
    Vec3f b1 = cross(xyz(q), b0);
    return p + 2 * (b0 * q.w + b1);
}

inline Quatf QuatConj(const Quatf& q)
{
    return Quatf(-q.x, -q.y, -q.z, q.w);
}

inline Quatf SLerp(const Quatf& q1, const Quatf& q2, float s)
{
    float cosHalfTheta = dot(q1, q2);

    if (abs(cosHalfTheta) >= float(0.99999))
        return q1;

    float sinHalfTheta = sqrt(float(1) - cosHalfTheta * cosHalfTheta);

    if (sinHalfTheta < float(1e-5))
        return float(0.5) * (q1 + q2);

    float halfTheta = atan2f(sinHalfTheta, cosHalfTheta);

    float t = float(1) - s;
    float ratio1 = sinf(t * halfTheta) / sinHalfTheta;
    float ratio2 = sinf(s * halfTheta) / sinHalfTheta;

    return ratio1 * q1 + ratio2 * q2;
}

inline Vec3f LogUnit3(const Quatf& q)
{
    float s = len(xyz(q));
    float c = q.w;
    return (atan2f(s, c) / (s + float(1e-8))) * xyz(q);
}

inline Quatf ExpUnit3(const Vec3f& lq)
{
    float theta = len(lq);
    Quatf q;
    xyz(q) = lq * (sinf(theta) / (theta + float(1e-8)));
    q.w    =       cosf(theta);
    return q;
}

inline Vec3f QuatDiff3(const Quatf& a, const Quatf& b)
{
    Vec3f v;
    v.x = + a.w * b.x - a.z * b.y + a.y * b.z - a.x * b.w;
    v.y = + a.z * b.x + a.w * b.y - a.x * b.z - a.y * b.w;
    v.z = - a.y * b.x + a.x * b.y + a.w * b.z - a.z * b.w;

    float c = dot(a, b);
    float s = len(v);

    v *= (atan2f(s, c) / (s + float(1e-8)));

    return v;
}

#endif
