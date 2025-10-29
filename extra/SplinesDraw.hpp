//
// SplinesDraw.hpp
//
// Draw support for splines
//
// Andrew Willmott
//

#ifndef SL_SPLINES_DRAW_H
#define SL_SPLINES_DRAW_H

#include "Colour.hpp"
#include "Splines.hpp"

namespace SL
{
    class Draw;

    enum SplineDrawMode
    {
        kSplineDrawModeNone,
        kSplineDrawModeT,           // Show t value
        kSplineDrawModePosition,    // Show position as rgb
        kSplineDrawModeVelocity,    // Show ||velocity|| as green->red
        kSplineDrawModeAcceleration,// Show ||acceleration|| as green->red
        kSplineDrawModeCurvature,   // Show curvature
        kSplineDrawModeLength,      // Show length from green (0) to red (1)
        kSplineDrawModeLine,        // Show lines that make up the rendered splines
        kNumSplineDrawModes
    };

    // 2D
    void DrawSplines          (Draw* dd, int numSplines, const Spline2 splines[], float tolerance = 0.05f);
    void DrawSplines          (Draw* dd, int numSplines, const Spline2 splines[], SplineDrawMode drawMode, float drawScale = 1.0f, float tolerance = 0.05f);

    void DrawSplinesWithWidth (Draw* dd, int numSplines, const Spline2 splines[], float width, float tolerance = 0.05f); // draw two curves representing splines with given width
    void FillSplinesWithWidth (Draw* dd, int numSplines, const Spline2 splines[], float width, float tolerance = 0.05f); // fill area of splines with given width

    void DrawSplinePoints     (Draw* dd, int numSplines, const Spline2 splines[], float r);
    void DrawSplineTangents   (Draw* dd, int numSplines, const Spline2 splines[], const Colour& colour = kColourBlue);
    void DrawSplineFrames     (Draw* dd, int numSplines, const Spline2 splines[], float size, float step = 0.2f, bool consistent = true);

    void DrawSplineCurvature  (Draw* dd, int numSplines, const Spline2 splines[], float step = 0.2f, float maxR = 100);

    void DrawSplinesMesh      (Draw* dd, int numSplines, const Spline2 splines[], float width, float tolerance = 0.05f);

    void DrawSplineExactBounds(Draw* dd, int numSplines, const Spline2 splines[], const Colour& colour = kColourOrange);
    void DrawSplineFastBounds (Draw* dd, int numSplines, const Spline2 splines[], const Colour& colour = kColourYellow);

    // 3D
    void DrawSplines          (Draw* dd, int numSplines, const Spline3 splines[], float tolerance = 0.05f);
    void DrawSplines          (Draw* dd, int numSplines, const Spline3 splines[], SplineDrawMode drawMode, float drawScale = 1.0f, float tolerance = 0.05f);

    void DrawSplinePoints     (Draw* dd, int numSplines, const Spline3 splines[], float r, Vec3f normal = vl_z);
    void DrawSplineTangents   (Draw* dd, int numSplines, const Spline3 splines[], const Colour& colour = kColourBlue);
    void DrawSplineFrames     (Draw* dd, int numSplines, const Spline3 splines[], float size, float step = 0.2f, bool consistent = true);

    void DrawSplineCurvature  (Draw* dd, int numSplines, const Spline3 splines[], float step = 0.2f, float maxR = 100);
    void DrawSplineTorsion    (Draw* dd, int numSplines, const Spline3 splines[], float step = 0.2f, float scale = 100);

    void DrawSplinesMesh      (Draw* dd, int numSplines, const Spline3 splines[], float width, float tolerance = 0.05f);

    void DrawSplinesWithWidth (Draw* dd, int numSplines, const Spline3 splines[], float width, float tolerance = 0.05f); // draw two curves representing splines with given width
    void FillSplinesWithWidth (Draw* dd, int numSplines, const Spline3 splines[], float width, float tolerance = 0.05f); // fill area of splines with given width

    void DrawSplineExactBounds(Draw* dd, int numSplines, const Spline3 splines[], const Colour& colour = kColourOrange);
    void DrawSplineFastBounds (Draw* dd, int numSplines, const Spline3 splines[], const Colour& colour = kColourYellow);

    // Helpers if you want to do your own spline->line conversions
    struct DLContext
    {
        Draw*          dd;
        SplineDrawMode drawMode;
        float          drawModeScale;
    };

    // Adapters for SplinesToLines*
    void DrawSplineLines           (int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* dlContext);
    void DrawSplineLinesWithColours(int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* dlContext);
    void DrawSplineLines           (int numLines, Vec3f lines[][2], const Spline3& spline, int i, int n, float ts[][2], void* dlContext);
    void DrawSplineLinesWithColours(int numLines, Vec3f lines[][2], const Spline3& spline, int i, int n, float ts[][2], void* dlContext);

    // Draw 2D splines using 3D api
    void DrawSplineLines3D           (int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* dlContext);
    void DrawSplineLinesWithColours3D(int numLines, Vec2f lines[][2], const Spline2& spline, int i, int n, float ts[][2], void* dlContext);
}

#endif
