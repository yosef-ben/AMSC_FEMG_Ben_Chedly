cmake_minimum_required(VERSION 3.16)
project(AMSC_FEMG_Ben_Chedly CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O2")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

set(BOOST_ROOT "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/boost/1.76.0")
set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost 1.76.0 REQUIRED COMPONENTS graph)

# ── Libreria FEMG ─────────────────────────────────────────────────────────────
add_library(femg
  FEMG/src/quantum_graph_problem.cpp
  FEMG/src/elliptic_problem.cpp
  FEMG/src/parabolic_problem.cpp
)

target_include_directories(femg PUBLIC
  ${CMAKE_SOURCE_DIR}/FEMG/include
)

target_link_libraries(femg PUBLIC
  Eigen3::Eigen
  Boost::graph
)

# ── Test parabolico ───────────────────────────────────────────────────────────
add_executable(test_heat
  test_problem/parabolic/test_heat.cpp
)

target_link_libraries(test_heat PRIVATE femg)

# ── Test ellittico ────────────────────────────────────────────────────────────
add_executable(test_elliptic
  test_problem/elliptic/test_elliptic.cpp
)

target_link_libraries(test_elliptic PRIVATE femg)
