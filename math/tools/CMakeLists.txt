define_module(BINARY math:tools
  DEPENDS dbglog>=1.7 utility>=1.47 service
  math>=${math_VERSION})

add_executable(math.io io.cpp)
buildsys_binary(math.io)
target_link_libraries(math.io ${MODULE_LIBRARIES})
target_compile_definitions(math.io PRIVATE ${MODULE_DEFINITIONS})
