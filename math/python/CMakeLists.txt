define_module(LIBRARY pymath=${MODULE_math_VERSION}
  DEPENDS
  math=${MODULE_math_VERSION}
  pysupport>=1.6
  Boost_THREAD Boost_SYSTEM)

set(pymath_SOURCES
  mathmodule.cpp
  )

add_library(pymath STATIC ${pymath_SOURCES})
buildsys_library(pymath)
target_link_libraries(pymath ${MODULE_LIBRARIES})
