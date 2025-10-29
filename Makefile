CXXFLAGS = --std=c++11

SPLINE_LIB_INCLUDES=Splines.hpp VLMini.hpp
SPLINE_LIB_SOURCES=Splines.cpp

splines_test: $(SPLINE_LIB_INCLUDES) $(SPLINE_LIB_SOURCES) SplinesTest.cpp
	$(CXX) $(CXXFLAGS) -o $@ $(SPLINE_LIB_SOURCES) SplinesTest.cpp

clean:
	$(RM) splines_test

test: splines_test
	@./splines_test > test.txt
	@! (diff test.txt test-ref.txt) || echo "All good"
