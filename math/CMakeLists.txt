define_module(LIBRARY math DEPENDS dbglog utility Boost)

set(math_SOURCES
	math.cpp math.hpp
    geometry_core.hpp
	geometry.cpp geometry.hpp
	filters.cpp filters.hpp
	combinatorics.cpp combinatorics.hpp
	signal.cpp signal.hpp
    transform.hpp
)

add_library(math STATIC ${math_SOURCES})
buildsys_library(math)
target_link_libraries(math ${MODULE_LIBRARIES} -lm)
