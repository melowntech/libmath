define_module(LIBRARY math DEPENDS dbglog Boost)

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
set_target_properties(math PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
  ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(math ${MODULE_LIBRARIES} -lm)
