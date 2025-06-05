####################### COMPARISON ############################

cmake_minimum_required(VERSION 3.3)

# set mode (Debug or Release)
set(CMAKE_BUILD_TYPE Release)

# set C++ standard
set(CMAKE_CXX_STANDARD 17)

project(MPC_MIQP_Solver_Comparison C)

option(CXX "enable C++ compilation" ON)

if(CXX)
  enable_language(CXX)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

find_package(GUROBI REQUIRED)
find_package(Threads REQUIRED)

include_directories(include 
    examples/maps
    examples/maps/utilities
    ../extern/polypartition/src
    ../mpc-classes/include
    ../mi-solver/include
    ../extern/eigen 
    ../qp-custom/cpp/include
    ../extern/zonocpp/include
    ${GUROBI_INCLUDE_DIRS}
)

add_executable(${CMAKE_PROJECT_NAME}
   examples/double_integrator_obs_avoid_solver_comparison.cpp
   src/MPC_MIQP_Hrep.cpp
   src/MPC_MIQP_Gurobi.cpp
   src/MPC_MIQP_Hrep_Gurobi.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_Comparison_3D
   examples/double_integrator_3D_obs_avoid_solver_comparison.cpp
   src/MPC_MIQP_Hrep.cpp
   src/MPC_MIQP_Gurobi.cpp
   src/MPC_MIQP_Hrep_Gurobi.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_Baseline
    examples/double_integrator_obs_avoid_baseline.cpp
    ../mpc-classes/src/MPC_MIQP.cpp
    ../mpc-classes/src/MPC_QP.cpp
    ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
    ../qp-custom/cpp/src/QP_primal_dual.cpp
    ../mi-solver/src/Data.cpp
    ../mi-solver/src/LeavesQueue.cpp
    ../mi-solver/src/Node.cpp
    ../mi-solver/src/MI_MPC.cpp
    ../mi-solver/src/DataVolatile.cpp
    ../mi-solver/src/VerboseOutputs.cpp
    ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_Baseline_3D
    examples/double_integrator_3D_obs_avoid_baseline.cpp
    ../mpc-classes/src/MPC_MIQP.cpp
    ../mpc-classes/src/MPC_QP.cpp
    ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
    ../qp-custom/cpp/src/QP_primal_dual.cpp
    ../mi-solver/src/Data.cpp
    ../mi-solver/src/LeavesQueue.cpp
    ../mi-solver/src/Node.cpp
    ../mi-solver/src/MI_MPC.cpp
    ../mi-solver/src/DataVolatile.cpp
    ../mi-solver/src/VerboseOutputs.cpp
    ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_NoDiag
   examples/double_integrator_obs_avoid_nodiag.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_NoDiag_3D
   examples/double_integrator_3D_obs_avoid_nodiag.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_NoReach
   examples/double_integrator_obs_avoid_noreach.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_NoReach_3D
   examples/double_integrator_3D_obs_avoid_noreach.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_NoDiag_NoReach
   examples/double_integrator_obs_avoid_nodiag_noreach.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_NoDiag_NoReach_3D
   examples/double_integrator_3D_obs_avoid_nodiag_noreach.cpp
   ../mpc-classes/src/MPC_MIQP.cpp
   ../mpc-classes/src/MPC_QP.cpp
   ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
   ../qp-custom/cpp/src/QP_primal_dual.cpp
   ../mi-solver/src/Data.cpp
   ../mi-solver/src/LeavesQueue.cpp
   ../mi-solver/src/Node.cpp
   ../mi-solver/src/MI_MPC.cpp
   ../mi-solver/src/DataVolatile.cpp
   ../mi-solver/src/VerboseOutputs.cpp
   ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_1_Thread
        examples/double_integrator_obs_avoid_1_thread.cpp
        ../mpc-classes/src/MPC_MIQP.cpp
        ../mpc-classes/src/MPC_QP.cpp
        ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
        ../qp-custom/cpp/src/QP_primal_dual.cpp
        ../mi-solver/src/Data.cpp
        ../mi-solver/src/LeavesQueue.cpp
        ../mi-solver/src/Node.cpp
        ../mi-solver/src/MI_MPC.cpp
        ../mi-solver/src/DataVolatile.cpp
        ../mi-solver/src/VerboseOutputs.cpp
        ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_8_Thread
        examples/double_integrator_obs_avoid_8_thread.cpp
        ../mpc-classes/src/MPC_MIQP.cpp
        ../mpc-classes/src/MPC_QP.cpp
        ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
        ../qp-custom/cpp/src/QP_primal_dual.cpp
        ../mi-solver/src/Data.cpp
        ../mi-solver/src/LeavesQueue.cpp
        ../mi-solver/src/Node.cpp
        ../mi-solver/src/MI_MPC.cpp
        ../mi-solver/src/DataVolatile.cpp
        ../mi-solver/src/VerboseOutputs.cpp
        ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_1_Thread_3D
        examples/double_integrator_3D_obs_avoid_1_thread.cpp
        ../mpc-classes/src/MPC_MIQP.cpp
        ../mpc-classes/src/MPC_QP.cpp
        ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
        ../qp-custom/cpp/src/QP_primal_dual.cpp
        ../mi-solver/src/Data.cpp
        ../mi-solver/src/LeavesQueue.cpp
        ../mi-solver/src/Node.cpp
        ../mi-solver/src/MI_MPC.cpp
        ../mi-solver/src/DataVolatile.cpp
        ../mi-solver/src/VerboseOutputs.cpp
        ../extern/polypartition/src/polypartition.cpp
)

add_executable(MPC_MIQP_Solver_8_Thread_3D
        examples/double_integrator_3D_obs_avoid_8_thread.cpp
        ../mpc-classes/src/MPC_MIQP.cpp
        ../mpc-classes/src/MPC_QP.cpp
        ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
        ../qp-custom/cpp/src/QP_primal_dual.cpp
        ../mi-solver/src/Data.cpp
        ../mi-solver/src/LeavesQueue.cpp
        ../mi-solver/src/Node.cpp
        ../mi-solver/src/MI_MPC.cpp
        ../mi-solver/src/DataVolatile.cpp
        ../mi-solver/src/VerboseOutputs.cpp
        ../extern/polypartition/src/polypartition.cpp
)

if(CXX)
 target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
                       debug ${GUROBI_CXX_DEBUG_LIBRARY})
 target_link_libraries(MPC_MIQP_Solver_Comparison_3D optimized ${GUROBI_CXX_LIBRARY}
                       debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_LIBRARY} Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_Comparison_3D ${GUROBI_LIBRARY} Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_Baseline Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_Baseline_3D Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_NoDiag Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_NoDiag_3D Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_NoReach Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_NoReach_3D Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_1_Thread Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_8_Thread Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_1_Thread_3D Threads::Threads)
target_link_libraries(MPC_MIQP_Solver_8_Thread_3D Threads::Threads)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()