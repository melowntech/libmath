# bump version here
set(math_VERSION 1.3)

define_module(LIBRARY math=${math_VERSION}
  DEPENDS dbglog>=1.4 utility>=1.11 Boost)

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
