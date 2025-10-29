#ifndef VL_MINI_H
#define VL_MINI_H

#include <math.h>

template<class T> inline T sqr(T x) { return x * x; }

template<class T> inline T vl_min(T a, T b) { return a < b ? a : b; }
template<class T> inline T vl_max(T a, T b) { return a > b ? a : b; }

constexpr float vlf_pi = 3.14159265358979323846f;
constexpr float vlf_halfPi = 0.5 * vlf_pi;
constexpr float vlf_twoPi = 2 * vlf_pi;

#define VI(T) T operator[](int i) const { return (&x)[i]; } T& operator[](int i){ return (&x)[i]; }

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

#endif
