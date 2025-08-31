
cmake_minimum_required(VERSION 3.16)
project(transient_heat)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(deal.II REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)

deal_ii_initialize_cached_variables()

add_executable(TransientThermal TransientThermal.cpp)
set_target_properties(TransientThermal PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
deal_ii_setup_target(TransientThermal)
