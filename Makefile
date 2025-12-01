CXXFLAGS = --std=c++11 -Wall -O2

SPLINE_LIB_INCLUDES=Splines.hpp RotSplines.hpp VLMini.hpp Makefile
SPLINE_LIB_SOURCES=Splines.cpp RotSplines.cpp

splines_test: $(SPLINE_LIB_INCLUDES) $(SPLINE_LIB_SOURCES) SplinesTest.cpp
	$(CXX) $(CXXFLAGS) -o $@ $(SPLINE_LIB_SOURCES) SplinesTest.cpp

clean:
	$(RM) splines_test

test: splines_test
	@./splines_test > test.txt
	@! (diff test.txt test-ref.txt) || echo "All good"
