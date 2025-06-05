# enable use of diaongal matrix in QP solver
cat ../qp-custom/cpp/include/diag_flag_0.hpp > ../qp-custom/cpp/include/diag_flag.hpp

# compile and run the simulation study
cd build
cmake ..
cmake --build .
./MPC_MIQP_Solver_Baseline
./MPC_MIQP_Solver_Baseline_3D
./MPC_MIQP_Solver_NoReach
./MPC_MIQP_Solver_NoReach_3D
./MPC_MIQP_Solver_1_Thread
./MPC_MIQP_Solver_8_Thread
./MPC_MIQP_Solver_1_Thread_3D
./MPC_MIQP_Solver_8_Thread_3D
cd ..

# disable use of diaongal matrix in QP solver
cat ../qp-custom/cpp/include/diag_flag_1.hpp > ../qp-custom/cpp/include/diag_flag.hpp

# compile and run the simulation study
cd build
cmake ..
cmake --build .
./MPC_MIQP_Solver_NoDiag
./MPC_MIQP_Solver_NoDiag_3D
./MPC_MIQP_Solver_NoDiag_NoReach
./MPC_MIQP_Solver_NoDiag_NoReach_3D
cd ..

# reset the diag_flag.hpp file
cat ../qp-custom/cpp/include/diag_flag_0.hpp > ../qp-custom/cpp/include/diag_flag.hpp