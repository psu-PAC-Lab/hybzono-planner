#include "MPC_QP.hpp"
#include "MpcMathUtilities.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>
#include <cmath>

// Ackermann nonlinear dynamics
Eigen::Vector<double,6> ackermannDyn(const Eigen::Vector<double,6> &x, const Eigen::Vector<double,2> &u,
  double wheelBase, double dT)
{
  // state: x = [xp; yp; th; psi; v; psidot]
  // input: u = [a; psiddot]

  // unpack
  double xp, yp, th, psi, v, psidot, a, psiddot;

  xp = x[0];
  yp = x[1];
  th = x[2];
  psi = x[3];
  v = x[4];
  psidot = x[5];

  a = u[0];
  psiddot = u[1];

  // Ackermann kinematics
  double xpdot, ypdot, thdot;
  xpdot = v*cos(th);
  ypdot = v*sin(th);
  thdot = (1/wheelBase)*v*tan(psi);

  // assemble derivative
  Eigen::Vector<double, 6> xdot;
  xdot << xpdot, ypdot, thdot, psidot, a, psiddot;

  // discrete time dynamics
  return x + xdot*dT;

}

// linearized dynamics
void ackermannLinearizedDynMatrices(const Eigen::Vector<double,6> &x, double wheelBase, 
  Eigen::Matrix<double,6,6> &Ac, Eigen::Matrix<double,6,3> &Bc)
{
  // unpack
  double xp, yp, th, psi, v, psidot;

  xp = x[0];
  yp = x[1];
  th = x[2];
  psi = x[3];
  v = x[4];
  psidot = x[5];

  // A matrix
  Eigen::Matrix<double, 6, 6> A;
  A << 0, 0, -v*sin(th), 0, cos(th), 0, 
       0, 0, v*cos(th), 0, sin(th), 0,
       0, 0, 0, (v/wheelBase)*pow(1/cos(psi),2), (1/wheelBase)*tan(psi), 0,
       0, 0, 0, 0, 0, 1,
       0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0;
  Ac = A;  

  // B matrix
  Eigen::Matrix<double, 6, 3> B;
  B << 0, 0, v*th*sin(th),
       0, 0, -v*th*cos(th),
       0, 0, -(1/wheelBase)*v*psi*pow(1/cos(psi), 2),
       0, 0, 0,
       1, 0, 0,
       0, 1, 0;
  Bc = B;

}

// function to construct a reference
std::vector<Eigen::VectorXd> constructReference(const Eigen::VectorXd &x0, int n_horizon, unsigned int n_sim, double dT)
{
  // constants
  const double pi = M_PI;

  // initialize output
  std::vector<Eigen::VectorXd> x_ref (n_horizon, Eigen::VectorXd::Zero(6));

  // rotational and velocity frequencies
  double om_th = 2*pi/10;
  double om_v = 2*pi/10;

  // oscillation mean and magnitude
  double v_avg = 5;
  double v_amp = 1;
  double th_avg = 0;
  double th_amp = pi/4;

  // integrate to current sim step
  Eigen::VectorXd x_n = x0; // init
  double t, th, v;
  for (int j=0; j<n_sim; j++)
  {
    // get v and th
    t = j*dT;
    th = x0[2] + th_amp*sin(om_th*t);
    v = x0[4] + v_amp*sin(om_v*t);

    // integrate
    x_n(0) = x_n(0) + v*cos(th)*dT;
    x_n(1) = x_n(1) + v*sin(th)*dT;
    x_n(2) = th;
    x_n(3) = 0;
    x_n(4) = v;
    x_n(5) = 0;
  }

  // integrate over horizon
  for (int k=0; k<n_horizon; k++)
  {

    // get v and th
    t = (n_sim+k)*dT;
    th = x0[2] + th_amp*sin(om_th*t);
    v = x0[4] + v_amp*sin(om_v*t);

    // integrate
    if (k==0)
    {
        x_ref[k](0) = x_n(0) + v*cos(th)*dT;
        x_ref[k](1) = x_n(1) + v*sin(th)*dT;
    }
    else
    {
        x_ref[k](0) = x_ref[k-1](0) + v*cos(th)*dT;
        x_ref[k](1) = x_ref[k-1](1) + v*sin(th)*dT;
    }
    x_ref[k](2) = th;
    x_ref[k](3) = 0;
    x_ref[k](4) = v;
    x_ref[k](5) = 0; 
  }

  return x_ref;
}

int main()
{
    // constants
    const double inf = std::numeric_limits<double>::max();
    const double pi = M_PI;

    // sim parameters
    const double dT = 0.1;
    const int n_horizon = 20;

    /* Ackermann vehicle (LTV) */

    // parameters
    const double wheelBase = 1;

    // state constraints (H-rep)
    Eigen::Matrix<double, 6, 6> I6 = Eigen::Matrix<double, 6, 6>::Identity();
    Eigen::Matrix<double, 12, 6> Ax_ineq;
    Eigen::Vector<double, 12> bx_ineq;
    Ax_ineq << I6, -1*I6; 
    bx_ineq << 1e4, 1e4, 1e4, pi/4, 10, pi,
        1e4, 1e4, 1e4, pi/4, 5, pi;

    // state constraint slack cost (H-rep)
    Eigen::Matrix<double, 12, 12> Qx_ineq_cost = 1e6*Eigen::Matrix<double, 12, 12>::Identity();
    Eigen::Vector<double, 12> sigma_max_x = 1e4*Eigen::Vector<double, 12>::Ones();

    // input constraints (H-rep)
    Eigen::Matrix<double, 3, 3> I3 = Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Matrix<double, 6, 3> Au_ineq;
    Eigen::Vector<double, 6> bu_ineq;
    Au_ineq << I3, -1*I3;
    bu_ineq << 10, 10*pi, 1, 10, 10*pi, -1;   

    // // state constraints as conzono
    // Eigen::Matrix<double, 3, 3> G_x = Eigen::Matrix<double, 3, 3>::Zero();
    // G_x.diagonal() << pi/4, 10, pi;
    // Eigen::Vector<double, 3> c_x = Eigen::Vector<double, 3>::Zero();
    // Eigen::Matrix<double, 0, 3> A_x = Eigen::Matrix<double, 0, 3>::Zero();
    // Eigen::Vector<double, 0> b_x = Eigen::Vector<double, 0>::Zero();
    // HybZonoMPC::ConZono Zc_x (G_x.sparseView(), c_x, A_x.sparseView(), b_x);
    // Eigen::Vector<double, 6> x_min_hard = -1e4*Eigen::Vector<double, 6>::Ones();
    // Eigen::Vector<double, 6> x_max_hard = 1e4*Eigen::Vector<double, 6>::Ones();
    // Eigen::Matrix<double, 3, 6> Pp_x;
    // Pp_x << 0, 0, 0, 1, 0, 0,
    //         0, 0, 0, 0, 1, 0,
    //         0, 0, 0, 0, 0, 1;

    // Eigen::Matrix<double, 3, 3> Qx_slack_cost = 1e6*Eigen::Matrix<double, 3, 3>::Identity();
    // Eigen::Vector<double, 3> sigma_max_x = 1e4*Eigen::Vector<double, 3>::Ones();

    // // input constraints as conzono
    // Eigen::Matrix<double, 2, 2> G_u = Eigen::Matrix<double, 2, 2>::Zero();
    // G_u.diagonal() << 10, 10*pi;
    // Eigen::Vector<double, 2> c_u = Eigen::Vector<double, 2>::Zero();
    // Eigen::Matrix<double, 0, 2> A_u = Eigen::Matrix<double, 0, 2>::Zero();
    // Eigen::Vector<double, 0> b_u = Eigen::Vector<double, 0>::Zero();
    // HybZonoMPC::ConZono Zc_u (G_u.sparseView(), c_u, A_u.sparseView(), b_u);
    // Eigen::Vector<double, 3> u_min_hard, u_max_hard;
    // u_min_hard << -1e4, -1e4, 1;
    // u_max_hard << 1e4, 1e4, 1;
    // Eigen::Matrix<double, 2, 3> Pp_u;
    // Pp_u << 1, 0, 0,
    //         0, 1, 0;

    // cost function matrices
    Eigen::SparseMatrix<double> Q (6,6);
    Eigen::SparseMatrix<double> R (3,3);
    Eigen::SparseMatrix<double> P (6,6);

    Q.diagonal() << 1, 1, 0, 0, 0, 0;
    R.diagonal() << 1e-3, 1e-3*pow(180/pi,2), 0;
    P.diagonal() << 1, 1, 1*pow(180/pi,2), 0, 0, 0;

    // initial state and input
    Eigen::Vector<double, 6> x0;
    x0 << 0, 0, 0, 0, 5, 0;
    Eigen::Vector<double, 3> u0;
    u0 << 0, 0, 1;

    // initial reference
    std::vector<Eigen::VectorXd> x_ref;
    x_ref = constructReference(x0, n_horizon, 0, dT);

    // initial dynamic linearization (about reference)
    std::vector<Eigen::SparseMatrix<double>> A_vec;
    std::vector<Eigen::SparseMatrix<double>> B_vec;
      
    Eigen::Matrix<double,6,6> Ac;
    Eigen::Matrix<double,6,3> Bc;
    Eigen::MatrixXd Ad, Bd; // need to be delcared as MatrixXd for cont2DiscDynMatrices

    for (int k=0; k<n_horizon; k++)
    {
      // linearize
      ackermannLinearizedDynMatrices(x_ref[k], wheelBase, Ac, Bc);

      // discretize
      cont2DiscDynMatrices(Ac, Bc, dT, Ad, Bd);
      A_vec.push_back(Ad.sparseView());
      B_vec.push_back(Bd.sparseView());
    }

    // create object
    HybZonoMPC::MPC_QP mpc;

    // define problem
    mpc.set_dyn_matrices(A_vec, B_vec);
    mpc.set_stage_cost(Q, R);
    mpc.set_terminal_cost(P);
    
    // set state and input inequality constraints (H-rep)
    mpc.set_state_inequality_constraints(Ax_ineq.sparseView(), bx_ineq);
    mpc.set_input_inequality_constraints(Au_ineq.sparseView(), bu_ineq);
    mpc.set_state_inequality_constraints_slack_vars_cost(Qx_ineq_cost.sparseView(), sigma_max_x);

    // // set state and input inequality constraints (conzono)
    // mpc.set_state_inequality_constraints(Zc_x, Pp_x.sparseView(), x_min_hard, x_max_hard);
    // mpc.set_input_inequality_constraints(Zc_u, Pp_u.sparseView(), u_min_hard, u_max_hard);
    // mpc.set_state_inequality_constraints_slack_vars_cost(Qx_slack_cost.sparseView(), sigma_max_x);
    
    // settings
    HybZonoMPC::MPC_Settings mpc_settings;
    mpc_settings.n_horizon = n_horizon;
    mpc_settings.u1_control = true;
    mpc.set_mpc_settings(mpc_settings);

    // settings
    QP_IP_MPC::QP_settings qp_settings;
    qp_settings.mu_init = 1e7;
    qp_settings.mu_term = 1e-3;
    qp_settings.mu_max = 1e10;
    qp_settings.T_max = dT;
    mpc.set_qp_settings(qp_settings);

    // build controller
    mpc.build_controller();

    // init control and state variables
    Eigen::VectorXd x = x0;
    Eigen::VectorXd u = u0;
    Eigen::VectorXd u1 = u0;

    // time execution
    double tot_time = 0;
    double max_time = 0;

    // log sim data
    std::ofstream refFile("../examples/data/ackermannRefData.txt");
    std::ofstream stateFile("../examples/data/ackermannStateData.txt");
    std::ofstream inputFile("../examples/data/ackermannInputData.txt");
    std::ofstream tFile("../examples/data/ackermannTimeData.txt");

    // run MPC controller in loop and print control inputs to screen
    const int n_cycles = 1000;
    for (int i = 0; i < n_cycles; i++)
    {
        // update control
        u = u1;

        // generate reference
        x_ref = constructReference(x0, n_horizon, i, dT);

        // linearize about reference
        for (int k=0; k<n_horizon; k++)
        {
            // linearize
            ackermannLinearizedDynMatrices(x_ref[k], wheelBase, Ac, Bc);

            // discretize
            cont2DiscDynMatrices(Ac, Bc, dT, Ad, Bd);
            A_vec[k] = Ad.sparseView();
            B_vec[k] = Bd.sparseView();
        }

        // start timer
        auto start_cycle = std::chrono::high_resolution_clock::now();

        // get control input
        auto [u1_sol, feasible] = mpc.control(x, x_ref, u, A_vec, B_vec);
        u1 = u1_sol;

        // log execution time in seconds
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_cycle);
        max_time = std::max(max_time, duration.count()/1e6);
        tot_time += duration.count()/1e6;

        // print control input
        //std::cout << "u = " << u.transpose() << std::endl;
        std::cout << "J = " << mpc.get_objective() << std::endl;

        // step dynamics
        x = ackermannDyn(x, u.segment(0,2), wheelBase, dT);

        // data logging
        tFile << i*dT << std::endl;;

        for (int j=0; j<6; j++)
            stateFile << x(j) << " ";
        stateFile << std::endl;

        for (int k=0; k<n_horizon; k++)
        {
            for (int j=0; j<6; j++)
                refFile << x_ref[k](j) << " ";
            refFile << std::endl;
        }
        refFile << std::endl;

        for (int j=0; j<2; j++)
            inputFile << u(j) << " ";
        inputFile << std::endl;
    }

    // print max execution time
    std::cout << "Max execution time: " << max_time << " sec" << std::endl;
    std::cout << "Average execution time: " << tot_time/n_cycles << " sec" << std::endl;

    // close files
    refFile.close();
    stateFile.close();
    inputFile.close();
    tFile.close();

    return 0;
}