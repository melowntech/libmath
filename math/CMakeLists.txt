# bump version here
set(math_VERSION 1.4)

set(math_EXTRA_DEPENDS)
set(math_EXTRA_DEFINITIONS)

if (OpenCV_FOUND)
  message(STATUS "math: compiling in OpenCV support")
  list(APPEND math_EXTRA_DEFINITIONS MATH_HAS_OPENCV=1)
  list(APPEND math_EXTRA_DEPENDS OpenCV)
else()
  message(STATUS "math: compiling without OpenCV support")
endif()


define_module(LIBRARY math=${math_VERSION}
  DEPENDS dbglog>=1.4 utility>=1.11 Boost
  ${math_EXTRA_DEPENDS}
  DEFINITIONS ${math_EXTRA_DEFINITIONS})

set(math_SOURCES
	math.cpp math.hpp
    geometry_core.hpp
	geometry.cpp geometry.hpp
	filters.cpp filters.hpp
	combinatorics.cpp combinatorics.hpp
	signal.cpp signal.hpp
    transform.hpp
    io.hpp
)

add_library(math STATIC ${math_SOURCES})
buildsys_library(math)
target_link_libraries(math ${MODULE_LIBRARIES} -lm)
buildsys_target_compile_definitions(math ${MODULE_DEFINITIONS})

