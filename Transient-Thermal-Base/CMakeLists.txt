
cmake_minimum_required(VERSION 3.16)
project(transient_heat)

find_package(deal.II REQUIRED)
deal_ii_initialize_cached_variables()

add_executable(SingleElement SingleElement.cpp)
deal_ii_setup_target(SingleElement)
