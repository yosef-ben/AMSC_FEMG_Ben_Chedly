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

# FEMG library
add_library(femg
  FEMG/src/quantum_graph_problem.cpp
  FEMG/src/parabolic_problem.cpp
  FEMG/src/fisher_kolmogorov_problem.cpp
  FEMG/src/spectral_problem.cpp
)

target_include_directories(femg PUBLIC
  ${CMAKE_SOURCE_DIR}/FEMG/include
)

target_link_libraries(femg PUBLIC
  Eigen3::Eigen
  Boost::graph
)

# Sequential performance instrumentation (ordering study)
add_executable(test_ordering_study
  test_problem/performance/test_ordering_study.cpp
)
target_link_libraries(test_ordering_study PRIVATE femg)

add_executable(test_condensation
  test_problem/performance/test_condensation.cpp
)
target_link_libraries(test_condensation PRIVATE femg)

# OpenMP applies to this executable alone; the library and every other
# target build exactly as before.
find_package(OpenMP COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
  add_executable(test_condensation_omp
    test_problem/performance/test_condensation_omp.cpp
  )
  target_link_libraries(test_condensation_omp PRIVATE femg OpenMP::OpenMP_CXX)
endif()

# MPI applies to this executable alone, like OpenMP above.
find_package(MPI COMPONENTS CXX)
if(MPI_CXX_FOUND)
  add_executable(test_condensation_mpi
    test_problem/performance/test_condensation_mpi.cpp
  )
  target_link_libraries(test_condensation_mpi PRIVATE femg MPI::MPI_CXX)
  if(OpenMP_CXX_FOUND)
    add_executable(test_condensation_hybrid
      test_problem/performance/test_condensation_mpi.cpp
    )
    target_link_libraries(test_condensation_hybrid
      PRIVATE femg MPI::MPI_CXX OpenMP::OpenMP_CXX)
  endif()
endif()

# Parabolic examples and validation tests
add_executable(test_heat
  test_problem/parabolic/test_heat.cpp
)
target_link_libraries(test_heat PRIVATE femg)

add_executable(test_heat_star_constant
  test_problem/parabolic/test_heat_star_constant.cpp
)
target_link_libraries(test_heat_star_constant PRIVATE femg)

add_executable(test_heat_star_linear
  test_problem/parabolic/test_heat_star_linear.cpp
)
target_link_libraries(test_heat_star_linear PRIVATE femg)

add_executable(test_heat_star_sine
  test_problem/parabolic/test_heat_star_sine.cpp
)
target_link_libraries(test_heat_star_sine PRIVATE femg)

add_executable(test_heat_star_sine_convergence
  test_problem/parabolic/test_heat_star_sine_convergence.cpp
)
target_link_libraries(test_heat_star_sine_convergence PRIVATE femg)

add_executable(test_heat_star_radial_decay
  test_problem/parabolic/test_heat_star_radial_decay.cpp
)
target_link_libraries(test_heat_star_radial_decay PRIVATE femg)

add_executable(test_heat_star_radial_decay_time_convergence
  test_problem/parabolic/test_heat_star_radial_decay_time_convergence.cpp
)
target_link_libraries(test_heat_star_radial_decay_time_convergence PRIVATE femg)

add_executable(test_heat_star_localized
  test_problem/parabolic/test_heat_star_localized.cpp
)
target_link_libraries(test_heat_star_localized PRIVATE femg)

add_executable(test_heat_graphene_eigenmode
  test_problem/parabolic/test_heat_graphene_eigenmode.cpp
)
target_link_libraries(test_heat_graphene_eigenmode PRIVATE femg)

add_executable(test_heat_tree_eigenmode
  test_problem/parabolic/test_heat_tree_eigenmode.cpp
)
target_link_libraries(test_heat_tree_eigenmode PRIVATE femg)

add_executable(test_heat_graphene_eigenmode_time_convergence
  test_problem/parabolic/test_heat_graphene_eigenmode_time_convergence.cpp
)
target_link_libraries(test_heat_graphene_eigenmode_time_convergence PRIVATE femg)

add_executable(test_heat_graphene_energy_decay
  test_problem/parabolic/test_heat_graphene_energy_decay.cpp
)
target_link_libraries(test_heat_graphene_energy_decay PRIVATE femg)

add_executable(test_fisher_kolmogorov_logistic
  test_problem/fisher_kolmogorov/test_fisher_kolmogorov_logistic.cpp
)
target_link_libraries(test_fisher_kolmogorov_logistic PRIVATE femg)

add_executable(test_fisher_kolmogorov_connectome
  test_problem/fisher_kolmogorov/test_fisher_kolmogorov_connectome.cpp
)
target_link_libraries(test_fisher_kolmogorov_connectome PRIVATE femg)

add_executable(test_fisher_kolmogorov_1d_sensitivity
  test_problem/fisher_kolmogorov/test_fisher_kolmogorov_1d_sensitivity.cpp
)
target_link_libraries(test_fisher_kolmogorov_1d_sensitivity PRIVATE femg)

add_executable(test_fisher_kolmogorov_fornari83
  test_problem/fisher_kolmogorov/test_fisher_kolmogorov_fornari83.cpp
)
target_link_libraries(test_fisher_kolmogorov_fornari83 PRIVATE femg)

add_executable(verify_fisher_kolmogorov
  test_problem/fisher_kolmogorov/verify_fisher_kolmogorov.cpp
)
target_link_libraries(verify_fisher_kolmogorov PRIVATE femg)

add_executable(test_fisher_kolmogorov_corti83
  test_problem/fisher_kolmogorov/test_fisher_kolmogorov_corti83.cpp
)
target_link_libraries(test_fisher_kolmogorov_corti83 PRIVATE femg)

add_executable(test_heat_graphene_complexity
  test_problem/parabolic/test_heat_graphene_complexity.cpp
)
target_link_libraries(test_heat_graphene_complexity PRIVATE femg)

add_executable(test_graph_setup
  test_problem/parabolic/test_graph_setup.cpp
)
target_link_libraries(test_graph_setup PRIVATE femg)

add_executable(test_spectral_star
  test_problem/spectral/test_spectral_star.cpp
)
target_link_libraries(test_spectral_star PRIVATE femg)

add_executable(test_spectral_graphene
  test_problem/spectral/test_spectral_graphene.cpp
)
target_link_libraries(test_spectral_graphene PRIVATE femg)

add_executable(test_spectral_tree
  test_problem/spectral/test_spectral_tree.cpp
)
target_link_libraries(test_spectral_tree PRIVATE femg)

add_executable(test_spectral_comparison
  test_problem/spectral/test_spectral_comparison.cpp
)
target_link_libraries(test_spectral_comparison PRIVATE femg)

add_executable(test_spectral_graphene_complexity
  test_problem/spectral/test_spectral_graphene_complexity.cpp
)
target_link_libraries(test_spectral_graphene_complexity PRIVATE femg)

add_executable(test_spectral_mesh_comparison
  test_problem/spectral/test_spectral_mesh_comparison.cpp
)
target_link_libraries(test_spectral_mesh_comparison PRIVATE femg)

add_executable(test_heat_convergence
  test_problem/parabolic/test_heat_convergence.cpp
)
target_link_libraries(test_heat_convergence PRIVATE femg)

add_executable(test_heat_time_convergence
  test_problem/parabolic/test_heat_time_convergence.cpp
)
target_link_libraries(test_heat_time_convergence PRIVATE femg)

add_executable(test_linear_solver_complexity
  test_problem/elliptic/test_linear_solver_complexity.cpp
)
target_link_libraries(test_linear_solver_complexity PRIVATE femg)

add_executable(test_fisher_kolmogorov_1d_time_step_study
  test_problem/fisher_kolmogorov/test_fisher_kolmogorov_1d_time_step_study.cpp
)
target_link_libraries(test_fisher_kolmogorov_1d_time_step_study PRIVATE femg)
