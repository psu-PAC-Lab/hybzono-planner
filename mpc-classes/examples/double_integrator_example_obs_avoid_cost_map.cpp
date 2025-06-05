#include "MPC_MIQP.hpp"
#include "ZonoCpp.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <iostream>
#include <chrono>
#include <map>
#include <vector>
#include <cmath>
#include <utility>
#include <fstream>

using namespace HybZonoMPC;
using namespace ZonoCpp;

void get_map_zonotope_tiling(std::pair<double,double> x_range, std::pair<double,double> y_range,
    HybZono &X, 
    std::map<int, Eigen::MatrixXd> &A_Hrep, std::map<int, Eigen::VectorXd> &b_Hrep,
    Eigen::VectorXd &cost_vec)
{
    // free space grid
    std::vector<std::vector<int>> free_space_grid;
    for (int i=0; i<12; i++)
    {
        free_space_grid.push_back({0,1,2,3,4});
    }

    std::vector<int> i_grid, j_grid;
    for (int i=0; i<free_space_grid.size(); i++)
    {
        for (int j : free_space_grid[i])
        {
            i_grid.push_back(i);
            j_grid.push_back(j);
        }
    }

    // grid dimensions
    int m_grid = free_space_grid.size();
    int n_grid = 5;
    double dx = (x_range.second - x_range.first)/m_grid;
    double dy = (y_range.second - y_range.first)/n_grid;

    // get tile centers
    std::vector<double> xc_tiles, yc_tiles;
    for (int i=0; i<i_grid.size(); i++)
    {
        xc_tiles.push_back(x_range.first + (((double) i_grid[i]) + 0.5)*dx);
        yc_tiles.push_back(y_range.first + (((double) j_grid[i]) + 0.5)*dy);
    }

    // get cost map
    double max_cost = 200;
    double sigma_gauss = 2*dx;
    cost_vec = Eigen::VectorXd::Zero(xc_tiles.size()); // init
    for (int i=0; i<xc_tiles.size(); i++)
    {
        double dist = std::sqrt(std::pow(xc_tiles[i], 2) + std::pow(yc_tiles[i], 2));
        cost_vec(i) = max_cost*std::exp(-std::pow(dist, 2)/(2*std::pow(sigma_gauss, 2)));
    }

    // directly construct hybzono

    // generator
    Eigen::Matrix<double, 2, 2> Gc = Eigen::Matrix<double, 2, 2>::Zero();
    Gc.diagonal() << dx, dy;
    Eigen::MatrixXd Gb = Eigen::MatrixXd::Zero(2, xc_tiles.size());
    for (int i=0; i<xc_tiles.size(); i++)
    {
        Gb(0, i) = xc_tiles[i];
        Gb(1, i) = yc_tiles[i];
    }
    Eigen::Vector<double, 2> c;
    c << -dx/2, -dy/2;

    // equality constraints
    Eigen::MatrixXd Ac = Eigen::MatrixXd::Zero(1, 2);
    Eigen::MatrixXd Ab = Eigen::MatrixXd::Ones(1, xc_tiles.size());
    Eigen::Vector<double, 1> b = Eigen::Vector<double, 1>::Ones();

    // set hybzono
    X.set(Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b, true);

    // get H-rep

    A_Hrep.clear();
    b_Hrep.clear();

    // polytope for template cell
    Eigen::Matrix<double, 4, 2> A_temp;
    A_temp << 1, 0,
             -1, 0,
              0, 1,
              0, -1;
    Eigen::Vector<double, 4> b_temp;
    b_temp << dx/2, dx/2, dy/2, dy/2;

    // get H-rep for each tile
    Eigen::Vector<double, 2> c_tile;
    for (int i=0; i<xc_tiles.size(); i++)
    {
        c_tile(0) = xc_tiles[i];
        c_tile(1) = yc_tiles[i];
        A_Hrep[i+1] = A_temp; // 1-based indexing
        b_Hrep[i+1] = b_temp + A_temp*c_tile; // 1-based indexing
    }
}

// reference
std::vector<Eigen::VectorXd> get_reference(int n_horizon)
{
    Eigen::VectorXd x_ref_i(4);
    x_ref_i << 15, 0, 0, 0;
    return std::vector<Eigen::VectorXd> (n_horizon, x_ref_i);
}

// print to file
std::ofstream file("../examples/data/VerboseOutputs.txt"); // global
void print_to_file(const std::string &str)
{
    file << str << std::endl;
}

int main()
{
    // problem definition
    const double max_speed = 1;
    const double dT = 1;
    
    Eigen::Matrix<double, 4, 4> A_dyn;
    Eigen::Matrix<double, 4, 2> B_dyn;
    A_dyn << 1, dT, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, dT,
             0, 0, 0, 1;
    B_dyn << 0, 0,
             dT, 0,
             0, 0,
             0, dT;

    // cost matrices
    Eigen::Matrix<double, 4, 4> Q_cost, P_cost;
    Eigen::Matrix<double, 2, 2> R_cost;
    Q_cost = Eigen::Matrix<double, 4, 4>::Zero();
    P_cost = Eigen::Matrix<double, 4, 4>::Zero();
    R_cost = Eigen::Matrix<double, 2, 2>::Zero();
    Q_cost.diagonal() << 0.1, 0, 0.1, 0;
    R_cost.diagonal() << 10, 10;
    P_cost.diagonal() << 10, 0, 10, 0;

    // state constraint matrices: x1 in [-1, 1], x2 in [-1, 1]
    Eigen::Matrix<double, 8, 4> Ax_ineq = Eigen::Matrix<double, 8, 4>::Zero();
    Eigen::Vector<double, 8> bx_ineq;
    Ax_ineq.block(0, 0, 4, 4) << Eigen::Matrix<double, 4, 4>::Identity();
    Ax_ineq.block(4, 0, 4, 4) << -Eigen::Matrix<double, 4, 4>::Identity();
    bx_ineq << 1e2, max_speed, 1e2, max_speed, 1e2, max_speed, 1e2, max_speed;

    // terminal state constraint matrices, same as state constraint matrices here
    Eigen::Matrix<double, 8, 4> Ax_term_ineq = Ax_ineq;
    Eigen::Vector<double, 8> bx_term_ineq = bx_ineq;

    // input constraint matrices: u in [-1, 1]
    Eigen::Matrix<double, 4, 2> Au_ineq = Eigen::Matrix<double, 4, 2>::Zero();
    Au_ineq.block(0, 0, 2, 2) << Eigen::Matrix<double, 2, 2>::Identity();
    Au_ineq.block(2, 0, 2, 2) << -Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Vector<double, 4> bu_ineq = Eigen::Vector<double, 4>::Ones();

    // constraint softening
    Eigen::Matrix<double, 8, 8> Q_slack_cost = 1e6*Eigen::Matrix<double, 8, 8>::Identity(8, 8); // slack cost
    Eigen::Vector<double, 8> sigma_max = 1e4*Eigen::Vector<double, 8>::Ones();

    // prediction horizon
    int n_horizon = 15;

    // loop time of controller
    double T_max = 1;

    // initial conditon
    Eigen::Vector<double, 4> x0;
    x0 << -15, 0, 0, 0;

    // create map
    HybZono X_hz;
    std::map<int, Eigen::MatrixXd> A_Hrep;
    std::map<int, Eigen::VectorXd> b_Hrep;
    std::pair<double,double> x_range = std::make_pair(-15, 15);
    std::pair<double,double> y_range = std::make_pair(-10, 10);
    Eigen::VectorXd cost_vec;
    get_map_zonotope_tiling(x_range, y_range, X_hz, A_Hrep, b_Hrep, cost_vec);
    
    // hybzono params
    Eigen::Matrix<double, 2, 2> Q_hz = 1e6*Eigen::Matrix<double, 2, 2>::Identity(); // slack cost
    Eigen::Matrix<double, 2, 4> Pc;
    Pc << 1, 0, 0, 0,
          0, 0, 1, 0;

    // create object
    MPC_MIQP mpc;

    // define problem
    mpc.set_dyn_matrices(A_dyn.sparseView(), B_dyn.sparseView());
    mpc.set_stage_cost(Q_cost.sparseView(), R_cost.sparseView());
    mpc.set_terminal_cost(P_cost.sparseView());
    mpc.set_state_inequality_constraints(Ax_ineq.sparseView(), bx_ineq);
    mpc.set_input_inequality_constraints(Au_ineq.sparseView(), bu_ineq);
    //mpc.set_terminal_state_inequality_constraints(Ax_term_ineq.sparseView(), bx_term_ineq);

    //mpc.set_state_inequality_constraints_slack_vars_cost(Q_slack_cost.sparseView(), sigma_max);

    mpc.setup_hybzono_Hrep_constraints(X_hz, A_Hrep, b_Hrep, Pc.sparseView(), Q_hz.sparseView(), -1e4, 1e4, cost_vec);
    
    mpc.set_distance_limit(1.001*(sqrt(2)*max_speed)*dT);

    // settings
    MI_QPIPMPC::MI_QPIPMPC_Settings mi_settings;
    mi_settings.n_threads = 1;
    mi_settings.T_max = 0; // infinite
    mi_settings.verbose = true;
    mpc.set_MI_settings(mi_settings);

    MPC_Settings mpc_settings;
    mpc_settings.warm_start = true;
    mpc_settings.n_horizon = n_horizon;
    mpc_settings.u1_control = true;
    mpc.set_mpc_settings(mpc_settings);

    // verbosity
    mpc.set_output_stream_fcn([](const std::string &str) { std::cout << str << std::endl; });
    //mpc.set_output_stream_fcn(&print_to_file);

    // u1 control
    Eigen::Vector<double, 2> u0 = Eigen::Vector<double, 2>::Zero();

    // build controller
    mpc.build_controller();

    // reference
    std::vector<Eigen::VectorXd> x_ref = get_reference(n_horizon);

    // simulate
    Eigen::VectorXd x = x0;
    Eigen::VectorXd u = u0;
    Eigen::VectorXd u1 = u0;

    // time execution
    double tot_time = 0;
    double max_time = 0;

    // run MPC controller in loop and print control inputs to screen
    const int n_cycles = 40;
    for (int i = 0; i < n_cycles; i++)
    {
        // update control
        u = u1;

        // start timer
        auto start_cycle = std::chrono::high_resolution_clock::now();

        // get control input
        auto [u1_sol, feasible] = mpc.control(x, x_ref, u);
        u1 = u1_sol;

        // log execution time in second
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start_cycle);
        max_time = std::max(max_time, duration.count()/1e6);
        tot_time += duration.count()/1e6;

        // update state
        x = A_dyn*x + B_dyn*u;

        // print control input
        std::cout << "u = " << u.transpose() << std::endl;
    }

    // get total execution

    // print max execution time
    std::cout << "Max execution time: " << max_time << " sec" << std::endl;
    std::cout << "Average execution time: " << tot_time/n_cycles << " sec" << std::endl;

    file.close();

    return 0;
}