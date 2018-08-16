SplineLib
=========

Library for manipulating 2D and 3D splines. Functionality included:

* Creation from Bezier, Hermite, and Catmull-Rom forms
* Creation of an array of splines from an array of points and tension value
* Evaluation of position, velocity, acceleration, and curvature
* Length measurement
* Finding bounds
* Subdivision
* Finding the closest point on a run of splines
* Helpers for advancing a point along a spline at some given velocity

To build and run the test app:

    c++ --std=c++11 Splines.cpp SplinesTest.cpp -o splines && ./splines


Examples
--------

Splines from Points:

![points](images/points.gif "Splines from Points")

Fast and Conservative Bounds:

![bounds](images/bounds.gif "Fast and Conservative Bounds")

Closest Point on Spline:

![pick](images/pick.gif "Closest Point on Spline")

Moving Points on Splines:

![agent](images/agent.gif "Moving Points on Splines")
