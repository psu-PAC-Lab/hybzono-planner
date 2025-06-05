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

// global vars
const double speed_lim = 1;

// map
void get_map_zonotope_tiling(std::pair<double,double> x_range, std::pair<double,double> y_range,
    HybZono &X, 
    std::map<int, Eigen::MatrixXd> &A_Hrep, std::map<int, Eigen::VectorXd> &b_Hrep,
    std::map<int, Eigen::MatrixXd> &A_Hrep_pos, std::map<int, Eigen::VectorXd> &b_Hrep_pos)
{
    // free space grid
    std::vector<std::vector<int>> free_space_grid;
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({1, 2});
    free_space_grid.push_back({1});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({1, 3});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});

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

    // directly construct hybzono

    // generator
    Eigen::Matrix<double, 4, 4> Gc = Eigen::Matrix<double, 4, 4>::Zero();
    Gc.diagonal() << dx, dy, 2*speed_lim, 2*speed_lim;
    Eigen::MatrixXd Gb = Eigen::MatrixXd::Zero(4, xc_tiles.size());
    for (int i=0; i<xc_tiles.size(); i++)
    {
        Gb(0, i) = xc_tiles[i];
        Gb(1, i) = yc_tiles[i];
    }
    Eigen::Vector<double, 4> c;
    c << -dx/2, -dy/2, -speed_lim, -speed_lim;

    // equality constraints
    Eigen::MatrixXd Ac = Eigen::MatrixXd::Zero(1, 4);
    Eigen::MatrixXd Ab = Eigen::MatrixXd::Ones(1, xc_tiles.size());
    Eigen::Vector<double, 1> b = Eigen::Vector<double, 1>::Ones();

    // set hybzono
    X.set(Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b, true);

    // get H-rep
    A_Hrep.clear();
    b_Hrep.clear();
    A_Hrep_pos.clear();
    b_Hrep_pos.clear();

    // polytope for template cell
    Eigen::Matrix<double, 8, 4> A_temp;
    A_temp << 1, 0, 0, 0,
             -1, 0, 0, 0,
              0, 1, 0, 0,
              0, -1, 0, 0,
              0, 0, 1, 0,
              0, 0, -1, 0,
              0, 0, 0, 1,
              0, 0, 0, -1;
    Eigen::Vector<double, 8> b_temp;
    b_temp << dx/2, dx/2, dy/2, dy/2, speed_lim, speed_lim, speed_lim, speed_lim;

    // get H-rep for each tile
    Eigen::Vector<double, 4> c_tile = Eigen::Vector<double, 4>::Zero();
    for (int i=0; i<xc_tiles.size(); i++)
    {
        c_tile(0) = xc_tiles[i];
        c_tile(1) = yc_tiles[i];
        A_Hrep[i+1] = A_temp; // 1-based indexing
        b_Hrep[i+1] = b_temp + A_temp*c_tile; // 1-based indexing
    }

    // position H-rep

    // polytope for template cell
    Eigen::Matrix<double, 4, 2> A_pos_temp;
    A_pos_temp << 1, 0,
             -1, 0,
              0, 1,
              0, -1;
    Eigen::Vector<double, 4> b_pos_temp;
    b_pos_temp << dx/2, dx/2, dy/2, dy/2;

    // get H-rep for each tile
    Eigen::Vector<double, 2> c_pos_tile = Eigen::Vector<double, 2>::Zero();
    for (int i=0; i<xc_tiles.size(); i++)
    {
        c_pos_tile(0) = xc_tiles[i];
        c_pos_tile(1) = yc_tiles[i];
        A_Hrep_pos[i+1] = A_pos_temp; // 1-based indexing
        b_Hrep_pos[i+1] = b_pos_temp + A_pos_temp*c_pos_tile; // 1-based indexing
    }
}

// reference
std::vector<Eigen::VectorXd> get_reference(int n_horizon)
{
    Eigen::VectorXd x_ref_i(4);
    x_ref_i << 15, 0, 0, 0;
    return std::vector<Eigen::VectorXd> (n_horizon, x_ref_i);
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
    R_cost.diagonal() << 10, 10;
    P_cost.diagonal() << 10, 0, 10, 0;

    // // H-rep constraints

    // // state constraint matrices: x1 in [-1, 1], x2 in [-1, 1]
    // Eigen::Matrix<double, 8, 4> Ax_ineq = Eigen::Matrix<double, 8, 4>::Zero();
    // Eigen::Vector<double, 8> bx_ineq;
    // Ax_ineq.block(0, 0, 4, 4) << Eigen::Matrix<double, 4, 4>::Identity();
    // Ax_ineq.block(4, 0, 4, 4) << -Eigen::Matrix<double, 4, 4>::Identity();
    // bx_ineq << 1e2, max_speed, 1e2, max_speed, 1e2, max_speed, 1e2, max_speed;

    // // terminal state constraint matrices, same as state constraint matrices here
    // Eigen::Matrix<double, 8, 4> Ax_term_ineq = Ax_ineq;
    // Eigen::Vector<double, 8> bx_term_ineq = bx_ineq;

    // // input constraint matrices: u in [-1, 1]
    // Eigen::Matrix<double, 4, 2> Au_ineq = Eigen::Matrix<double, 4, 2>::Zero();
    // Au_ineq.block(0, 0, 2, 2) << Eigen::Matrix<double, 2, 2>::Identity();
    // Au_ineq.block(2, 0, 2, 2) << -Eigen::Matrix<double, 2, 2>::Identity();
    // Eigen::Vector<double, 4> bu_ineq = Eigen::Vector<double, 4>::Ones();

    // conzono constraints
    
    // state constraints - only on velocity
    Eigen::Matrix<double, 2, 2> Gx = Eigen::Matrix<double, 2, 2>::Zero();
    Gx.diagonal() << max_speed, max_speed;
    Eigen::Vector<double, 2> cx = Eigen::Vector<double, 2>::Zero();
    Eigen::Matrix<double, 0, 2> Ax = Eigen::Matrix<double, 0, 2>::Zero();
    Eigen::Vector<double, 0> bx = Eigen::Vector<double, 0>::Zero();
    ConZono Zc_x (Gx.sparseView(), cx, Ax.sparseView(), bx);
    Eigen::Vector<double, 4> x_max_hard = 1e4*Eigen::Vector<double, 4>::Ones();
    Eigen::Vector<double, 4> x_min_hard = -1e4*Eigen::Vector<double, 4>::Ones();
    Eigen::Matrix<double, 2, 4> Pp_x;
    Pp_x << 0, 1, 0, 0,
            0, 0, 0, 1;
    
    // soft state constraints
    Eigen::Matrix<double, 2, 2> Qx_ineq_cost = 1e6*Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Vector<double, 2> sigma_max_x = 1e4*Eigen::Vector<double, 2>::Ones();

    // input constraints
    Eigen::Matrix<double, 2, 2> Gu;
    Gu.diagonal() << 1, 1;
    Eigen::Vector<double, 2> cu = Eigen::Vector<double, 2>::Zero();
    Eigen::Matrix<double, 0, 2> Au = Eigen::Matrix<double, 0, 2>::Zero();
    Eigen::Vector<double, 0> bu = Eigen::Vector<double, 0>::Zero();
    ConZono Zc_u (Gu.sparseView(), cu, Au.sparseView(), bu);
    Eigen::Vector<double, 2> u_max_hard = 1e4*Eigen::Vector<double, 2>::Ones();
    Eigen::Vector<double, 2> u_min_hard = -1e4*Eigen::Vector<double, 2>::Ones();
    Eigen::Matrix<double, 2, 2> Pp_u = Eigen::Matrix<double, 2, 2>::Identity();


    // prediction horizon
    int n_horizon = 15;

    // loop time of controller
    double T_max = 1;

    // initial conditon
    Eigen::Vector<double, 4> x0;
    x0 << -15, 0, 0, 0;

    // create map
    HybZono X_hz;
    std::map<int, Eigen::MatrixXd> A_Hrep, A_Hrep_pos;
    std::map<int, Eigen::VectorXd> b_Hrep, b_Hrep_pos;
    std::pair<double,double> x_range = std::make_pair(-15, 15);
    std::pair<double,double> y_range = std::make_pair(-10, 10);
    get_map_zonotope_tiling(x_range, y_range, X_hz, A_Hrep, b_Hrep, A_Hrep_pos, b_Hrep_pos);
    
    // hybzono params
    Eigen::Matrix<double, 2, 4> Pp;
    Pp << 1, 0, 0, 0,
          0, 0, 1, 0;

    Eigen::Matrix<double, 4, 4> Pc;
    Pc << 1, 0, 0, 0,
          0, 0, 1, 0,
          0, 1, 0, 0,
          0, 0, 0, 1;

    Eigen::Matrix<double, 4, 4> Q_hz = 1e6*Eigen::Matrix<double, 4, 4>::Identity(); // slack cost

    // create object
    MPC_MIQP mpc;

    // define problem
    
    mpc.set_dyn_matrices(A_dyn.sparseView(), B_dyn.sparseView());
    mpc.set_stage_cost(Q_cost.sparseView(), R_cost.sparseView());
    mpc.set_terminal_cost(P_cost.sparseView());

    // // H-rep constraints
    // mpc.set_state_inequality_constraints(Ax_ineq.sparseView(), bx_ineq);
    // mpc.set_input_inequality_constraints(Au_ineq.sparseView(), bu_ineq);
    // mpc.set_terminal_state_inequality_constraints(Ax_term_ineq.sparseView(), bx_term_ineq);

    // conzono constraints
    mpc.set_state_inequality_constraints(Zc_x, Pp_x.sparseView(), x_min_hard, x_max_hard);
    mpc.set_input_inequality_constraints(Zc_u, Pp_u.sparseView(), u_min_hard, u_max_hard);
    //mpc.set_terminal_state_inequality_constraints(Zc_x, Pp_x.sparseView(), x_min_hard, x_max_hard); 
    mpc.set_state_inequality_constraints_slack_vars_cost(Qx_ineq_cost.sparseView(), sigma_max_x);

    mpc.setup_hybzono_Hrep_constraints(X_hz, A_Hrep, b_Hrep, Pc.sparseView(), Q_hz.sparseView(), -1e4, 1e4, Eigen::VectorXd(), 
        Pp.sparseView(), A_Hrep_pos, b_Hrep_pos);
    
    mpc.set_distance_limit(1.001*(sqrt(2)*max_speed)*dT);

    // settings
    MI_QPIPMPC::MI_QPIPMPC_Settings mi_settings;
    mi_settings.n_threads = 6;
    mi_settings.T_max = 0; // infinite
    mi_settings.verbose = false;
    mi_settings.conv_abs = 0;
    mi_settings.conv_rel = 0;
    mpc.set_MI_settings(mi_settings);

    MPC_Settings mpc_settings;
    mpc_settings.warm_start = false;
    mpc_settings.n_horizon = n_horizon;
    mpc_settings.u1_control = true;
    mpc.set_mpc_settings(mpc_settings);

    // verbosity
    mpc.set_output_stream_fcn([](const std::string &str) { std::cout << str << std::endl; });

    // u1 control
    Eigen::Vector<double, 2> u0 = Eigen::Vector<double, 2>::Zero();

    // build controller
    mpc.build_controller();

    // reference
    std::vector<Eigen::VectorXd> x_ref = get_reference(n_horizon);

    // output to file
    std::ofstream file;
    file.open("../examples/data/double_integrator_obs_avoid_example.txt");
    file << x0(0) << " " << x0(2) << std::endl;

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

        // log data
        file << x(0) << " " << x(2) << std::endl;
    }

    // get total execution

    // print max execution time
    std::cout << "Max execution time: " << max_time << " sec" << std::endl;
    std::cout << "Average execution time: " << tot_time/n_cycles << " sec" << std::endl;

    file.close();

    return 0;
}