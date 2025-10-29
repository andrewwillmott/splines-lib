SplineLib
=========

Library for manipulating 1D, 2D, and 3D splines. Functionality included:

* Creation from Bezier, Hermite, and Catmull-Rom forms
* Creation of an array of splines from an array of points and tension value,
  or Bezier hull points, or Hermite tangents.
* Creation of arcs and circles
* Evaluation of position, velocity, acceleration, curvature, torsion, and
  tangent frame
* Length measurement
* Finding bounds
* Offsetting (e.g., for stroking) and reversing splines
* Splitting and re-joining of single splines
* Subdivision of runs of splines either evenly, for flatness, or to be linear in
  arc length
* Conversion of runs of splines to line segments, e.g., for drawing
* Finding the closest point on a run of splines
* Finding where runs of splines intersect, or a run of splines self-intersects
* Helpers for advancing a point along a spline at some given velocity
* Support for monotonic splines, which don't extrapolate outside the source
  points

To build and run the test app:

    make test

Or, add the Splines* files to your favourite IDE.

In `extra` there is some code for drawing various spline features (as per the
examples below). It's missing the draw and transform dependencies, but could be
useful as reference.


Examples
--------

Splines from Points:

![points](images/points.gif "Splines from Points")

Fast and Conservative Bounds:

![bounds](images/bounds.gif "Fast and Conservative Bounds")

Closest Point on Spline:

![pick](images/pick.gif "Closest Point on Spline")

Spline Intersections:

![pick](images/self-intersect.gif "Spline Intersections")

Moving Points on Splines:

![agent](images/agent.gif "Moving Points on Splines")

Subdivision and Joining:

![pick](images/subdivide.gif "Subdivision and Joining")
