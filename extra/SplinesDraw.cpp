//
// SplinesDraw.cpp
//
// Draw support for splines
//
// Andrew Willmott
//

#include "SplinesDraw.hpp"

#include "Draw.hpp"
#include "MinMax.hpp"
#include "Transform2.hpp"
#include "Transform3.hpp"

using namespace SL;

namespace
{
    inline float Fract(float x)
    {
        return x - floorf(x);
    }

    Colour Fract(const Colour& c)
    {
        return { Fract(c[0]), Fract(c[1]), Fract(c[2]) };
    }

    Colour SplineColour(const Spline2& spline, float t, SplineDrawMode mode, float scale, int i)
    {
        switch (mode)
        {
        case kSplineDrawModeNone:
            scale = 0.0f;
            break;

        case kSplineDrawModeT:
            scale *= t;
            break;

        case kSplineDrawModePosition:
            return Fract(scale * Vec3f(Position(spline, t), 0.0));

        case kSplineDrawModeVelocity:
            scale *= len(Velocity(spline, t)) * 0.001f;
            break;

        case kSplineDrawModeAcceleration:
            scale *= len(Acceleration(spline, t)) * 0.01f;

        case kSplineDrawModeCurvature:
            scale *= Curvature(spline, t) * 100.0f;
            break;

        case kSplineDrawModeLength:
            if (t >= 0.5f)
                scale *= Length(spline, t, 1.0f) * 0.01f;
            else
                scale *= Length(spline, 0.0f, t) * 0.01f;
            break;

        case kSplineDrawModeLine:
            return LabelColour(i);
            break;

        case kNumSplineDrawModes:
            SL_ERROR("bad enum");
        }

        scale = ClampUnit(scale);
        float hue = lerp(kHueGreen, kHueRed + 360.0f, scale);   // move through the blues to red

        return HSVColour(hue, 1.0f, 1.0f);
    }
}

void SL::DrawSplineLines(int numLines, Vec2f lines[][2], const Spline2&, int i, int n, float[][2], void* context)
{
    reinterpret_cast<Draw*>(context)->DrawLines(2 * numLines, lines[0]);
}

void SL::DrawSplineLinesWithColours(int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* contextIn)
{
    DLContext& ctx = *reinterpret_cast<DLContext*>(contextIn);

    ColourAlpha c[kMaxEmitLines][2];

    for (int i = 0; i < numLines; i++)
    {
        c[i][0] = SplineColour(spline, ts[i][0], ctx.drawMode, ctx.drawModeScale, i);
        c[i][1] = SplineColour(spline, ts[i][1], ctx.drawMode, ctx.drawModeScale, i);
    }

    ctx.dd->DrawLines(numLines * 2, lines[0], c[0]);
}

void SL::DrawSplines(Draw* dd, int numSplines, const Spline2 splines[], float tolerance)
{
    SplinesToLinesAdaptive(numSplines, splines, ::DrawSplineLines, dd, tolerance);
}

void SL::DrawSplines(Draw* dd, int numSplines, const Spline2 splines[], SplineDrawMode drawMode, float drawScale, float tolerance)
{
    DLContext ctx = { dd, drawMode, drawScale };
    SplinesToLinesAdaptive(numSplines, splines, ::DrawSplineLinesWithColours, &ctx, tolerance);
}

namespace
{
    void OffsetLR(const Spline2& spline, float offset, Spline2* spline_0, Spline2* spline_1)
    {
        offset *= 0.5f;

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

        Vec4f ox(sy0 * sd0, sy1 * sd1, sy2 * sd2, sy3 * sd3);
        Vec4f oy(sx0 * sd0, sx1 * sd1, sx2 * sd2, sx3 * sd3);

        *spline_0 = {spline.xb - ox, spline.yb + oy};
        *spline_1 = {spline.xb + ox, spline.yb - oy};
    }

    void OffsetLR(const Spline3& spline, float offset, Spline3* spline_0, Spline3* spline_1)
    {
        OffsetLR((const Spline2&) spline, offset, (Spline2*) spline_0, (Spline2*) spline_1);
        spline_0->zb = spline.zb;
        spline_1->zb = spline.zb;
    }

    struct DLContextWidth
    {
        Draw* dd;
        float width;
    };

    void DrawLinesWidth(int numLines, Vec3f lines[][2], const Spline3& spline, int, int, float ts[][2], void* contextIn)
    {
        DLContextWidth& context = *reinterpret_cast<DLContextWidth*>(contextIn);

        Vec3f p[kMaxEmitLines][4];

        Spline3 os0, os1;
        OffsetLR(spline, context.width, &os0, &os1);

        for (int i = 0; i < numLines; i++)
        {
            p[i][0] = Position(os0, ts[i][0]);
            p[i][1] = Position(os0, ts[i][1]);
            p[i][2] = Position(os1, ts[i][0]);
            p[i][3] = Position(os1, ts[i][1]);
        }

        context.dd->DrawLines(4 * numLines, p[0]);
    }

    void FillLinesWidth(int numLines, Vec3f lines[][2], const Spline3& spline, int, int, float ts[][2], void* contextIn)
    {
        DLContextWidth& context = *reinterpret_cast<DLContextWidth*>(contextIn);

        Vec3f p[kMaxEmitLines][12];

        Spline3 os0, os1;
        OffsetLR(spline, context.width, &os0, &os1);

        for (int i = 0; i < numLines; i++)
        {
            Vec3f v00 = Position(os0, ts[i][0]);
            Vec3f v01 = Position(os1, ts[i][0]);
            Vec3f v10 = Position(os0, ts[i][1]);
            Vec3f v11 = Position(os1, ts[i][1]);

            p[i][0] = v00;
            p[i][1] = v01;
            p[i][2] = v11;
            p[i][3] = v00;
            p[i][4] = v11;
            p[i][5] = v10;

            v00 = cross(norm_safe(lines[i][0]-v00), norm_safe(Velocity(spline, ts[i][0])));
            v10 = cross(norm_safe(lines[i][1]-v10), norm_safe(Velocity(spline, ts[i][1])));

            v01 = lines[i][0] - v00 * context.width * 0.5f;
            v00 = lines[i][0] + v00 * context.width * 0.5f;
            v11 = lines[i][1] - v10 * context.width * 0.5f;
            v10 = lines[i][1] + v10 * context.width * 0.5f;

            p[i][6] = v00;
            p[i][7] = v01;
            p[i][8] = v11;
            p[i][9] = v00;
            p[i][10] = v11;
            p[i][11] = v10;
        }

        context.dd->DrawTriangles(12 * numLines, p[0]);
    }

    void DrawLinesWidth(int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* contextIn)
    {
        DLContextWidth& context = *reinterpret_cast<DLContextWidth*>(contextIn);

        Vec2f p[kMaxEmitLines][4];

        Spline2 os0, os1;
        OffsetLR(spline, context.width, &os0, &os1);

        for (int i = 0; i < numLines; i++)
        {
            p[i][0] = Position(os0, ts[i][0]);
            p[i][1] = Position(os0, ts[i][1]);
            p[i][2] = Position(os1, ts[i][0]);
            p[i][3] = Position(os1, ts[i][1]);
        }

        context.dd->DrawLines(4 * numLines, p[0]);
    }

    void FillLinesWidth(int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* contextIn)
    {
        DLContextWidth& context = *reinterpret_cast<DLContextWidth*>(contextIn);

        Vec2f p[kMaxEmitLines][6];

        Spline2 os0, os1;
        OffsetLR(spline, context.width, &os0, &os1);

        for (int i = 0; i < numLines; i++)
        {
            Vec2f v00 = Position(os0, ts[i][0]);
            Vec2f v01 = Position(os1, ts[i][0]);
            Vec2f v10 = Position(os0, ts[i][1]);
            Vec2f v11 = Position(os1, ts[i][1]);

            p[i][0] = v00;
            p[i][1] = v01;
            p[i][2] = v11;
            p[i][3] = v00;
            p[i][4] = v11;
            p[i][5] = v10;
        }

        context.dd->DrawTriangles(6 * numLines, p[0]);
    }
}


void SL::DrawSplinesWithWidth(Draw* dd, int numSplines, const Spline2 splines[], float width, float tolerance)
{
    DLContextWidth context = {dd, width};
    SplinesToLinesAdaptive(numSplines, splines, ::DrawLinesWidth, &context, tolerance);
}

void SL::FillSplinesWithWidth(Draw* dd, int numSplines, const Spline2 splines[], float width, float tolerance)
{
    DLContextWidth context = {dd, width};
    SplinesToLinesAdaptive(numSplines, splines, ::FillLinesWidth, &context, tolerance);
}

void SL::DrawSplinesWithWidth(Draw* dd, int numSplines, const Spline3 splines[], float width, float tolerance)
{
    DLContextWidth context = {dd, width};
    SplinesToLinesAdaptive(numSplines, splines, ::DrawLinesWidth, &context, tolerance);
}

void SL::FillSplinesWithWidth(Draw* dd, int numSplines, const Spline3 splines[], float width, float tolerance)
{
    DLContextWidth context = {dd, width};
    SplinesToLinesAdaptive(numSplines, splines, ::FillLinesWidth, &context, tolerance);
}

//////////////////////////////////////////////////


void SL::DrawSplinePoints(Draw* dd, int numSplines, const Spline2 splines[], float r)
{
    if (numSplines > 0)
    {
        dd->SetColour(kColourPurple);

        FillCircle(dd, Position0(splines[0             ]), r);
        FillCircle(dd, Position1(splines[numSplines - 1]), r);
    }

    r *= 0.5f;
    dd->SetColour(kColourRed);

    for (int i = 1; i < numSplines; i++)
        FillCircle(dd, Position0(splines[i]), r);
}

void SL::DrawSplineTangents(Draw* dd, int numSplines, const Spline2 splines[], const Colour& colour)
{
    dd->SetColour(colour);

    for (int i = 0; i < numSplines; i++)
    {
        const Spline2& spline = splines[i];

        Vec2f p0 = Position0(spline);
        Vec2f p1 = Position1(spline);

        Vec2f v0 = Velocity0(spline);
        Vec2f v1 = Velocity1(spline);

        Vec2f p[] =
        {
            p0,
            p0 + v0 / 3.0f,
            p1,
            p1 - v1 / 3.0f,
        };

        dd->DrawLines(SL_SIZE(p), p);
    }
}

void SL::DrawSplineFrames(Draw* dd, int numSplines, const Spline2 splines[], float size, float step, bool consistent)
{
    Mat2f F = vl_I;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline2& spline = splines[i];

        for (float t = 0.5f * step; t <= 1.0f; t += step)
        {
            Vec2f p = Position(spline, t);

            if (consistent)
                Frame(spline, t, &F);
            else
                F= Frame(spline, t);

            Transform2 xform(size, F, p);
            dd->PushTransform2D(xform);

            dd->SetColour(kColourRed);
            DrawLine(dd, Vec2f(vl_0), vl_x);
            dd->SetColour(kColourGreen);
            DrawLine(dd, Vec2f(vl_0), vl_y);

            dd->PopTransform2D();
        }
    }
}

void SL::DrawSplineCurvature(Draw* dd, int numSplines, const Spline2 splines[], float step, float maxR)
{
    float minK = 1.0f / maxR;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline2& spline = splines[i];

        for (float t = 0.0f; t <= 1.0f; t += step)
        {
            float k = Curvature(spline, t);

            if (k < minK)
                continue;

            Vec2f p = Position(spline, t);
            float r = 1.0f / k;
            Mat2f F = Frame(spline, t);

            DrawCircle(dd, p + F.y * r, r);
        }
    }
}

void SL::DrawSplineExactBounds(Draw* dd, int numSplines, const Spline2 splines[], const Colour& colour)
{
    dd->SetColour(colour);

    for (int i = 0; i < numSplines; i++)
        DrawRect(dd, ExactBounds(splines[i]));
}

void SL::DrawSplineFastBounds(Draw* dd, int numSplines, const Spline2 splines[], const Colour& colour)
{
    dd->SetColour(colour);

    for (int i = 0; i < numSplines; i++)
        DrawRect(dd, FastBounds(splines[i]));
}

// 3D

namespace
{
    Colour SplineColour(const Spline3& spline, float t, SplineDrawMode mode, float scale, int i)
    {
        switch (mode)
        {
        case kSplineDrawModeNone:
            scale = 0.0f;
            break;

        case kSplineDrawModeT:
            scale *= t;
            break;

        case kSplineDrawModePosition:
            return Fract(scale * Position(spline, t));

        case kSplineDrawModeVelocity:
            scale *= len(Velocity(spline, t)) * 0.001f;
            break;

        case kSplineDrawModeAcceleration:
            scale *= len(Acceleration(spline, t)) * 0.01f;

        case kSplineDrawModeCurvature:
            scale *= Curvature(spline, t) * 100.0f;
            break;

        case kSplineDrawModeLength:
            scale *= Length(spline, 0.0f, t) * 0.01f;
            break;

        case kSplineDrawModeLine:
            return LabelColour(i);
            break;

        case kNumSplineDrawModes:
            SL_ERROR("bad enum");
        }

        scale = ClampUnit(scale);
        float hue = lerp(kHueGreen, kHueRed + 360.0f, scale);

        return HSVColour(hue, 1.0f, 1.0f);
    }
}

void SL::DrawSplineLines(int numLines, Vec3f lines[][2], const Spline3&, int, int, float[][2], void* context)
{
    reinterpret_cast<Draw*>(context)->DrawLines(2 * numLines, lines[0]);
}

void SL::DrawSplineLinesWithColours(int numLines, Vec3f lines[][2], const Spline3& spline, int, int, float ts[][2], void* contextIn)
{
    DLContext& ctx = *reinterpret_cast<DLContext*>(contextIn);

    ColourAlpha c[kMaxEmitLines][2];

    for (int i = 0; i < numLines; i++)
    {
        c[i][0] = SplineColour(spline, ts[i][0], ctx.drawMode, ctx.drawModeScale, i);
        c[i][1] = SplineColour(spline, ts[i][1], ctx.drawMode, ctx.drawModeScale, i);
    }

    ctx.dd->DrawLines(numLines * 2, lines[0], c[0]);
}

void SL::DrawSplineLines3D(int numLines, Vec2f lines[][2], const Spline2&, int i, int n, float[][2], void* context)
{
    Draw* dd = reinterpret_cast<Draw*>(context);

    Vec3f p3[kMaxEmitLines][2];

    for (int i = 0; i < numLines; i++)
    {
        p3[i][0] = { lines[i][0].x, lines[i][0].y, 0.01f };
        p3[i][1] = { lines[i][1].x, lines[i][1].y, 0.01f };
    }

    dd->DrawLines(2 * numLines, p3[0]);
}

void SL::DrawSplineLinesWithColours3D(int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* contextIn)
{
    DLContext& ctx = *reinterpret_cast<DLContext*>(contextIn);

    Vec3f p3[kMaxEmitLines][2];
    ColourAlpha c[kMaxEmitLines][2];

    for (int i = 0; i < numLines; i++)
    {
        p3[i][0] = { lines[i][0].x, lines[i][0].y, 0.01f };
        p3[i][1] = { lines[i][1].x, lines[i][1].y, 0.01f };

        c[i][0] = SplineColour(spline, ts[i][0], ctx.drawMode, ctx.drawModeScale, i);
        c[i][1] = SplineColour(spline, ts[i][1], ctx.drawMode, ctx.drawModeScale, i);
    }

    ctx.dd->DrawLines(numLines * 2, p3[0], c[0]);
}

void SL::DrawSplines(Draw* dd, int numSplines, const Spline3 splines[], float tolerance)
{
    SplinesToLinesAdaptive(numSplines, splines, ::DrawSplineLines, dd, tolerance);
}

void SL::DrawSplines(Draw* dd, int numSplines, const Spline3 splines[], SplineDrawMode drawMode, float drawScale, float tolerance)
{
    DLContext ctx = { dd, SplineDrawMode(drawMode), drawScale };
    SplinesToLinesAdaptive(numSplines, splines, ::DrawSplineLinesWithColours, &ctx, tolerance);
}

void SL::DrawSplinePoints(Draw* dd, int numSplines, const Spline3 splines[], float r, Vec3f normal)
{
    if (numSplines > 0)
    {
        dd->SetColour(kColourPurple);

        FillCircle(dd, Position0(splines[0             ]), r, normal);
        FillCircle(dd, Position1(splines[numSplines - 1]), r, normal);
    }

    r *= 0.5f;
    dd->SetColour(kColourRed);

    for (int i = 1; i < numSplines; i++)
        FillCircle(dd, Position0(splines[i]), r, normal);
}

void SL::DrawSplineTangents(Draw* dd, int numSplines, const Spline3 splines[], const Colour& colour)
{
    dd->SetColour(colour);

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        Vec3f p0 = Position0(spline);
        Vec3f p1 = Position1(spline);

        Vec3f v0 = Velocity0(spline);
        Vec3f v1 = Velocity1(spline);

        Vec3f p[] =
        {
            p0,
            p0 + v0 / 3.0f,
            p1,
            p1 - v1 / 3.0f,
        };

        dd->DrawLines(SL_SIZE(p), p);
    }
}

void SL::DrawSplineFrames(Draw* dd, int numSplines, const Spline3 splines[], float size, float step, bool consistent)
{
    Mat3f F = vl_I;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        for (float t = 0.5f * step; t <= 1.0f; t += step)
        {
            Vec3f p = Position(spline, t);

            if (consistent)
                Frame(spline, t, &F);
            else
                F = Frame(spline, t);

            Transform3 xform(size, F, p);
            dd->PushTransform3D(xform);

            dd->SetColour(kColourRed);
            DrawLine(dd, Vec3f(vl_0), vl_x);
            dd->SetColour(kColourGreen);
            DrawLine(dd, Vec3f(vl_0), vl_y);
            dd->SetColour(kColourBlue);
            DrawLine(dd, Vec3f(vl_0), vl_z);

            dd->PopTransform3D();
        }
    }
}

void SL::DrawSplineCurvature(Draw* dd, int numSplines, const Spline3 splines[], float step, float maxR)
{
    float minK = 1.0f / maxR;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        for (float t = 0.5f * step; t <= 1.0f; t += step)
        {
            float k = Curvature(spline, t);

            if (k < minK)
                continue;

            Vec3f p = Position(spline, t);
            float r = 1.0f / k;
            Mat3f F = Frame(spline, t);

            Transform3 xform(r, F, p);
            dd->PushTransform3D(xform);
            DrawCircle(dd, Vec3f(vl_y), 1.0f, Vec3f(vl_z));
            dd->PopTransform3D();
        }
    }
}

void SL::DrawSplineTorsion(Draw* dd, int numSplines, const Spline3 splines[], float step, float scale)
{
    Vec3f c0 = dd->Colour();
    Vec3f c1 = Vec3f(vl_1) - c0;

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        for (float t = 0.5f * step; t <= 1.0f; t += step)
        {
            float tt = Torsion(spline, t);
            float r = fabsf(tt * scale);

            if (r < 1e-1f)
                continue;

            Vec3f p = Position(spline, t);
            Mat3f F = Frame(spline, t);

            Transform3 xform(r, F, p);
            dd->PushTransform3D(xform);
            if (tt < 0.0f)
                dd->SetColour(c1);
            else
                dd->SetColour(c0);
            DrawCircle(dd, Vec3f(vl_0), 1.0f, Vec3f(vl_x));
            dd->PopTransform3D();
        }
    }

    dd->SetColour(c0);
}

void SL::DrawSplineExactBounds(Draw* dd, int numSplines, const Spline3 splines[], const Colour& colour)
{
    dd->SetColour(colour);

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        Bounds3 bounds = ExactBounds(spline);

        DrawBox(dd, bounds.mMin, bounds.mMax);
    }
}

void SL::DrawSplineFastBounds(Draw* dd, int numSplines, const Spline3 splines[], const Colour& colour)
{
    dd->SetColour(colour);

    for (int i = 0; i < numSplines; i++)
    {
        const Spline3& spline = splines[i];

        Bounds3 bounds = FastBounds(spline);

        DrawBox(dd, bounds.mMin, bounds.mMax);
    }
}
