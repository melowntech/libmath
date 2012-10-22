define_module(LIBRARY math DEPENDS dbglog Boost)

set(math_SOURCES
	math.cpp
	geometry.cpp
	filters.cpp
	combinatorics.cpp
	signal.cpp
)

add_library(math STATIC ${math_SOURCES})
set_target_properties(math PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
  ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(math ${MODULE_LIBRARIES} -lm)
