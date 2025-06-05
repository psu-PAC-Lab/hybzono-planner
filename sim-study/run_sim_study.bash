# enable use of diaongal matrix in QP solver
cat ../qp-custom/cpp/include/diag_flag_0.hpp > ../qp-custom/cpp/include/diag_flag.hpp

# compile and run the simulation study
cd build
cmake ..
cmake --build .
./MPC_MIQP_Solver_Comparison
./MPC_MIQP_Solver_Comparison_3D
cd ..