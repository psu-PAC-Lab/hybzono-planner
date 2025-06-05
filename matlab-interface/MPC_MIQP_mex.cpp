// WORK IN PROGRESS!
//
// class_wrapper_template.cpp
// Example of using a C++ class via a MEX-file
// by Jonathan Chappelow (chappjc)
//
// Design goals:
//   1. Manage multiple persistent instances of a C++ class
//   2. Small consecutive integer handles used in MATLAB (not cast pointers)
//   3. Transparently handle resource management (i.e. MATLAB never 
//      responsible for memory allocated for C++ classes)
//       a. No memory leaked if MATLAB fails to issue "delete" action
//       b. Automatic deallocation if MEX-file prematurely unloaded
//   4. Guard against premature module unloading
//   5. Validity of handles implicitly verified without checking a magic number
//   6. No wrapper class or functions mimicking mexFunction, just an intuitive
//      switch-case block in mexFunction.
//
// Note that these goals should be acheved without regard to any MATLAB class, 
// but which can also help address memory management issues.  As such, the
// resulting MEX-file can safely be used directly (but not too elegantly).
//
// Use:
//   1. Enumerate the different actions (e.g. New, Delete, Insert, etc.) in the
//      Actions enum.  For each enumerated action, specify a string (e.g.
//      "new", "delete", "insert", etc.) to be passed as the first argument to
//      the MEX function in MATLAB.
//   2. Customize the handling for each action in the switch statement in the
//      body of mexFunction (e.g. call the relevant C++ class method).
//  
// Implementation:
//
// For your C++ class, class_type, mexFunction uses static data storage to hold
// a persistent (between calls to mexFunction) table of integer handles and 
// smart pointers to dynamically allocated class instances.  A std::map is used
// for this purpose, which facilitates locating known handles, for which only 
// valid instances of your class are guaranteed to exist:
//
//    typedef unsigned int handle_type;
//    std::map<handle_type, std::shared_ptr<class_type>>
//
// A std::shared_ptr takes care of deallocation when either (1) a table element
// is erased via the "delete" action or (2) the MEX-file is unloaded.
//
// To prevent the MEX-file from unloading while a MATLAB class instances exist,
// mexLock is called each time a new C++ class instance is created, adding to
// the MEX-file's lock count.  Each time a C++ instance is deleted mexUnlock is
// called, removing one lock from the lock count.
//
// Requirements:
//
// A modern compiler with the following C++11 features:
//   - shared_ptr
//   - auto
//   - enum class
//   - initializer_list (for const map initialization)
// (VS2013, recent GCC possibly with -std=c++11, Clang since 3.1)
//
// TODO:
//
// - This example uses a priority queue class of mine, which is far from the
//   simplest example.  Demonstrate with something more basic.
// - Somehow put the relevant parts in a header, OR the other way around -- put
//   user class config in a header and include into the mexFunction's cpp.
//
// 04/25/15 (chappjc) - Initial version.

#include "mex.h"

#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <string>
#include <sstream>
#include <exception>


// JAR additions
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <stdlib.h>
#include "mex_set_methods.hpp"

using namespace MEXSET;

////////////////////////  BEGIN Step 1: Configuration  ////////////////////////

// Include your class declaration
#include "MPC_MIQP.hpp"

// Define class_type for your class
typedef double data_type; // DO I NEED THIS?
typedef HybZonoMPC::MPC_MIQP class_type;

// List actions
enum class Action
{
    // create/destroy instance - REQUIRED
    New,
    Delete,
    // user-specified class functionality
    Set_dyn_matrices,
    Set_dyn_matrices_LTV,
    Set_stage_cost,
    Set_terminal_cost,
    Set_state_inequality_constraints,
    Set_input_inequality_constraints,
    Set_terminal_state_inequality_constraints,
    Set_state_inequality_constraints_conzono,
    Set_input_inequality_constraints_conzono,
    Set_terminal_state_inequality_constraints_conzono,
    Set_state_inequality_constraints_slack_vars_cost,
    Set_input_inequality_constraints_slack_vars_cost,
    Set_terminal_state_inequality_constraints_slack_vars_cost,
    Set_QP_settings,
    Set_distance_limit,
    Setup_hybzono_Hrep_constraints,
    Setup_hybzono_constraints,
    Set_dist_between_regions,
    Set_terminal_hybzono_constraints,
    Set_MPC_settings,
    Set_MI_settings,
    Build_controller,
    Control,
    Control_LTV,
    Get_trajectory,
    Get_terminal_region_selection,
    Get_miqp_results,
    Set_region_dependent_disturbances
};

// Map string (first input argument to mexFunction) to an Action
const std::map<std::string, Action> actionTypeMap =
{
    { "new", Action::New },
    { "delete", Action::Delete },
    { "set_dyn_matrices", Action::Set_dyn_matrices },
    { "set_dyn_matrices_LTV", Action::Set_dyn_matrices_LTV },
    { "set_stage_cost", Action::Set_stage_cost },
    { "set_terminal_cost", Action::Set_terminal_cost },
    { "set_state_inequality_constraints", Action::Set_state_inequality_constraints },
    { "set_input_inequality_constraints", Action::Set_input_inequality_constraints },
    { "set_terminal_state_inequality_constraints", Action::Set_terminal_state_inequality_constraints },
    { "set_state_inequality_constraints_conzono", Action::Set_state_inequality_constraints_conzono },
    { "set_input_inequality_constraints_conzono", Action::Set_input_inequality_constraints_conzono },
    { "set_terminal_state_inequality_constraints_conzono", Action::Set_terminal_state_inequality_constraints_conzono },
    { "set_state_inequality_constraints_slack_vars_cost", Action::Set_state_inequality_constraints_slack_vars_cost },
    { "set_input_inequality_constraints_slack_vars_cost", Action::Set_input_inequality_constraints_slack_vars_cost },
    { "set_terminal_state_inequality_constraints_slack_vars_cost", Action::Set_terminal_state_inequality_constraints_slack_vars_cost },
    { "set_qp_settings", Action::Set_QP_settings },
    { "set_distance_limit", Action::Set_distance_limit },
    { "setup_hybzono_Hrep_constraints", Action::Setup_hybzono_Hrep_constraints },
    { "setup_hybzono_constraints", Action::Setup_hybzono_constraints },
    { "set_dist_between_regions", Action::Set_dist_between_regions },
    { "set_terminal_hybzono_constraints", Action::Set_terminal_hybzono_constraints },
    { "set_mpc_settings", Action::Set_MPC_settings },
    { "set_mi_settings", Action::Set_MI_settings },
    { "build_controller", Action::Build_controller },
    { "control", Action::Control },
    { "control_LTV", Action::Control_LTV },
    { "get_trajectory", Action::Get_trajectory },
    { "get_terminal_region_selection", Action::Get_terminal_region_selection },
    { "get_miqp_results", Action::Get_miqp_results },
    { "set_region_dependent_disturbances", Action::Set_region_dependent_disturbances }
};

// results struct field names
const char * MI_RESULTS_FIELDS[] = 
{ 
    "upper_glob",
    "lower_glob",
    "run_time",
    "status",
    "qp_solve_time",
    "qp_iter_avg",
    "iter_num"
};

// status key
double mi_status_2_double(MI_QPIPMPC::MI_Status status)
{
    switch (status)
    {
    case MI_QPIPMPC::MI_Status::MI_UNSOLVED:
        return 0;
    case MI_QPIPMPC::MI_Status::MI_SOLVED:
        return 1;
    case MI_QPIPMPC::MI_Status::MI_INFEASIBLE:
        return 2;
    case MI_QPIPMPC::MI_Status::MI_MAX_ITER_FEASIBLE:
        return 3;
    case MI_QPIPMPC::MI_Status::MI_MAX_ITER_UNSOLVED:
        return 4;
    default:
        return -1;
    }
}

/////////////////////////  END Step 1: Configuration  /////////////////////////

typedef unsigned int handle_type;
typedef std::pair<handle_type, std::shared_ptr<class_type>> indPtrPair_type; // or boost::shared_ptr
typedef std::map<indPtrPair_type::first_type, indPtrPair_type::second_type> instanceMap_type;
typedef indPtrPair_type::second_type instPtr_t;

// getHandle pulls the integer handle out of prhs[1]
handle_type getHandle(int nrhs, const mxArray *prhs[]);
// checkHandle gets the position in the instance table
instanceMap_type::const_iterator checkHandle(const instanceMap_type&, handle_type);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    // static storage duration object for table mapping handles to instances
    static instanceMap_type instanceTab;

    if (nrhs < 1 || !mxIsChar(prhs[0]))
        mexErrMsgTxt("First input must be an action string ('new', 'delete', or a method name).");

    char *actionCstr = mxArrayToString(prhs[0]); // convert char16_t to char
    std::string actionStr(actionCstr); mxFree(actionCstr);

    //for (auto & c : actionStr) c = ::tolower(c); // remove this for case sensitivity

    if (actionTypeMap.count(actionStr) == 0)
        mexErrMsgTxt(("Unrecognized action (not in actionTypeMap): " + actionStr).c_str());

    // If action is not "new" or "delete" try to locate an existing instance based on input handle
    instPtr_t instance;
    if (actionTypeMap.at(actionStr) != Action::New && actionTypeMap.at(actionStr) != Action::Delete) {
        handle_type h = getHandle(nrhs, prhs);
        instanceMap_type::const_iterator instIt = checkHandle(instanceTab, h);
        instance = instIt->second;
    }

    // verbosity
    MI_QPIPMPC::MI_QPIPMPC_Settings mi_settings_default;
    static bool verbose = mi_settings_default.verbose;
    std::string verbose_output; // declare

	//////// Step 2: customize the each action in the switch in mexFuction ////////
    switch (actionTypeMap.at(actionStr))
    {
    case Action::New:
    {
        if (nrhs > 1 && mxGetNumberOfElements(prhs[1]) != 1)
            mexErrMsgTxt("Second argument (optional) must be a scalar, N.");

        handle_type newHandle = instanceTab.size() ? (instanceTab.rbegin())->first + 1 : 1;

        std::pair<instanceMap_type::iterator, bool> insResult;
        insResult = instanceTab.insert(indPtrPair_type(newHandle, std::make_shared<class_type>())); // default constructor

        if (!insResult.second) // sanity check
            mexPrintf("Oh, bad news.  Tried to add an existing handle."); // shouldn't ever happen
        else
            mexLock(); // add to the lock count

		// return the handle
        plhs[0] = mxCreateDoubleScalar(insResult.first->first); // == newHandle

        break;
    }
    case Action::Delete:
    {
        instanceMap_type::const_iterator instIt = checkHandle(instanceTab, getHandle(nrhs, prhs));
        instanceTab.erase(instIt);
        mexUnlock();
        plhs[0] = mxCreateLogicalScalar(instanceTab.empty()); // info

        break;
    }
    case Action::Set_dyn_matrices:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two inputs required: handle, A_dyn, B_dyn.");

        // get input data pointers
        const mxArray *A_dyn_mx = prhs[2];
        const mxArray *B_dyn_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> A_dyn, B_dyn;
        std::vector<Eigen::Triplet<double>> A_dyn_triplets, B_dyn_triplets;
        A_dyn_triplets.reserve(getNNZ(A_dyn_mx));
        B_dyn_triplets.reserve(getNNZ(B_dyn_mx));
        getSparseMatrixInput(A_dyn_mx, A_dyn_triplets, &A_dyn);
        getSparseMatrixInput(B_dyn_mx, B_dyn_triplets, &B_dyn);

        // call method
        instance->set_dyn_matrices(A_dyn, B_dyn);

        break;
    }
    case Action::Set_dyn_matrices_LTV:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two inputs required: handle, A_dyn_vec, B_dyn_vec.");

        // get input data pointers
        const mxArray *A_dyn_vec_mx = prhs[2];
        const mxArray *B_dyn_vec_mx = prhs[3];

        // init input fields
        std::vector<Eigen::SparseMatrix<double>> A_dyn_vec, B_dyn_vec;
        std::vector<Eigen::Triplet<double>> A_dyn_vec_triplets, B_dyn_vec_triplets;
        Eigen::SparseMatrix<double> A_dyn_el, B_dyn_el;

        // loop through cell array to get matrices
        for (int i = 0; i < mxGetNumberOfElements(A_dyn_vec_mx); i++)
        {
            A_dyn_vec_triplets.clear();
            B_dyn_vec_triplets.clear();
            const mxArray *A_dyn_mx = mxGetCell(A_dyn_vec_mx, i);
            const mxArray *B_dyn_mx = mxGetCell(B_dyn_vec_mx, i);
            getSparseMatrixInput(A_dyn_mx, A_dyn_vec_triplets, &A_dyn_el);
            getSparseMatrixInput(B_dyn_mx, B_dyn_vec_triplets, &B_dyn_el);
            A_dyn_vec.push_back(A_dyn_el);
            B_dyn_vec.push_back(B_dyn_el);
        }

        // call method
        instance->set_dyn_matrices(A_dyn_vec, B_dyn_vec);

        break;
    }
    case Action::Set_stage_cost:
    {
        // check number of inputs
        if (nrhs != 5+2)
            mexErrMsgTxt("Five inputs required: Q, R, N, q_x, q_u.");

        // get input data pointers
        const mxArray *Q_mx = prhs[2];
        const mxArray *R_mx = prhs[3];
        const mxArray *N_mx = prhs[4];
        const mxArray *q_x_mx = prhs[5];
        const mxArray *q_u_mx = prhs[6];

        // get data
        Eigen::SparseMatrix<double> Q, R, N;
        std::vector<Eigen::Triplet<double>> Q_triplets, R_triplets, N_triplets;
        Q_triplets.reserve(getNNZ(Q_mx));
        R_triplets.reserve(getNNZ(R_mx));
        Eigen::VectorXd q_x, q_u;
        
        getSparseMatrixInput(Q_mx, Q_triplets, &Q);
        getSparseMatrixInput(R_mx, Q_triplets, &R);

        if (!mxIsEmpty(N_mx))
            getSparseMatrixInput(N_mx, N_triplets, &N);

        if (!mxIsEmpty(q_x_mx))
            getVectorInput(q_x_mx, q_x);

        if (!mxIsEmpty(q_x_mx))
            getVectorInput(q_u_mx, q_u);

        // call method
        instance->set_stage_cost(Q, R, N, q_x, q_u);

        break;
    }
    case Action::Set_terminal_cost:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two inputs required: P, q_xN");

        // get input data pointers
        const mxArray *P_mx = prhs[2];
        const mxArray *q_xN_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> P;
        std::vector<Eigen::Triplet<double>> P_triplets;
        P_triplets.reserve(getNNZ(P_mx));
        Eigen::VectorXd q_xN;

        getSparseMatrixInput(P_mx, P_triplets, &P);

        if (!mxIsEmpty(q_xN_mx))
            getVectorInput(q_xN_mx, q_xN);

        // call method
        instance->set_terminal_cost(P, q_xN);

        break;
    }
    case Action::Set_state_inequality_constraints:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two input required: handle, Ax_ineq, bx_ineq.");

        // get input data pointers
        const mxArray *Ax_ineq_mx = prhs[2];
        const mxArray *bx_ineq_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> Ax_ineq;
        Eigen::VectorXd bx_ineq;
        std::vector<Eigen::Triplet<double>> Ax_ineq_triplets;
        Ax_ineq_triplets.reserve(getNNZ(Ax_ineq_mx));
        getSparseMatrixInput(Ax_ineq_mx, Ax_ineq_triplets, &Ax_ineq);
        getVectorInput(bx_ineq_mx, bx_ineq);

        // call method
        instance->set_state_inequality_constraints(Ax_ineq, bx_ineq);

        break;
    }
    case Action::Set_input_inequality_constraints:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two input required: handle, Au_ineq, bu_ineq.");

        // get input data pointers
        const mxArray *Au_ineq_mx = prhs[2];
        const mxArray *bu_ineq_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> Au_ineq;
        Eigen::VectorXd bu_ineq;
        std::vector<Eigen::Triplet<double>> Au_ineq_triplets;
        Au_ineq_triplets.reserve(getNNZ(Au_ineq_mx));
        getSparseMatrixInput(Au_ineq_mx, Au_ineq_triplets, &Au_ineq);
        getVectorInput(bu_ineq_mx, bu_ineq);

        // call method
        instance->set_input_inequality_constraints(Au_ineq, bu_ineq);

        break;
    }
    case Action::Set_terminal_state_inequality_constraints:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two input required: handle, Ax_term_ineq, bx_term_ineq.");

        // get input data pointers
        const mxArray *Ax_term_ineq_mx = prhs[2];
        const mxArray *bx_term_ineq_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> Ax_term_ineq;
        Eigen::VectorXd bx_term_ineq;
        std::vector<Eigen::Triplet<double>> Ax_term_ineq_triplets;
        Ax_term_ineq_triplets.reserve(getNNZ(Ax_term_ineq_mx));
        getSparseMatrixInput(Ax_term_ineq_mx, Ax_term_ineq_triplets, &Ax_term_ineq);
        getVectorInput(bx_term_ineq_mx, bx_term_ineq);

        // call method
        instance->set_terminal_state_inequality_constraints(Ax_term_ineq, bx_term_ineq);

        break;
    }
    case Action::Set_state_inequality_constraints_conzono:
    {
        // check number of inputs
        if (nrhs != 2+4)
            mexErrMsgTxt("Four inputs required: handle, Zc_x, Pp_x, x_min_hard, x_max_hard.");

        // get input data pointers
        const mxArray *Zc_x_mx = prhs[2];
        const mxArray *Pp_x_mx = prhs[3];
        const mxArray *x_min_hard_mx = prhs[4];
        const mxArray *x_max_hard_mx = prhs[5];

        // declare
        ZonoCpp::ConZono Zc_x;
        Eigen::MatrixXd G_x, A_x, Pp_x;
        Eigen::VectorXd c_x, b_x, x_min_hard, x_max_hard;

        // get conzono
        getFullMatrixInput(mxGetField(Zc_x_mx, 0, "G"), G_x);
        getFullMatrixInput(mxGetField(Zc_x_mx, 0, "A"), A_x);
        getVectorInput(mxGetField(Zc_x_mx, 0, "c"), c_x);
        getVectorInput(mxGetField(Zc_x_mx, 0, "b"), b_x);
        Zc_x.set(G_x.sparseView(), c_x, A_x.sparseView(), b_x);

        // get Pp_x
        getFullMatrixInput(Pp_x_mx, Pp_x);

        // get x_min_hard and x_max_hard
        getVectorInput(x_min_hard_mx, x_min_hard);
        getVectorInput(x_max_hard_mx, x_max_hard);

        // call method
        instance->set_state_inequality_constraints(Zc_x, Pp_x.sparseView(), x_min_hard, x_max_hard);

        break;
    }
    case Action::Set_input_inequality_constraints_conzono:
    {
        // check number of inputs
        if (nrhs != 2+4)
            mexErrMsgTxt("Four inputs required: handle, Zc_u, Pp_u, u_min_hard, u_max_hard.");

        // get input data pointers
        const mxArray *Zc_u_mx = prhs[2];
        const mxArray *Pp_u_mx = prhs[3];
        const mxArray *u_min_hard_mx = prhs[4];
        const mxArray *u_max_hard_mx = prhs[5];

        // declare
        ZonoCpp::ConZono Zc_u;
        Eigen::MatrixXd G_u, A_u, Pp_u;
        Eigen::VectorXd c_u, b_u, u_min_hard, u_max_hard;

        // get conzono
        getFullMatrixInput(mxGetField(Zc_u_mx, 0, "G"), G_u);
        getFullMatrixInput(mxGetField(Zc_u_mx, 0, "A"), A_u);
        getVectorInput(mxGetField(Zc_u_mx, 0, "c"), c_u);
        getVectorInput(mxGetField(Zc_u_mx, 0, "b"), b_u);
        Zc_u.set(G_u.sparseView(), c_u, A_u.sparseView(), b_u);

        // get Pp_u
        getFullMatrixInput(Pp_u_mx, Pp_u);

        // get u_min_hard and u_max_hard
        getVectorInput(u_min_hard_mx, u_min_hard);
        getVectorInput(u_max_hard_mx, u_max_hard);

        // call method
        instance->set_input_inequality_constraints(Zc_u, Pp_u.sparseView(), u_min_hard, u_max_hard);

        break;
    }
    case Action::Set_terminal_state_inequality_constraints_conzono:
    {
        // check number of inputs
        if (nrhs != 2+4)
            mexErrMsgTxt("Four inputs required: handle, Zc_x_term, Pp_x_term, x_min_hard_term, x_max_hard_term.");

        // get input data pointers
        const mxArray *Zc_x_term_mx = prhs[2];
        const mxArray *Pp_x_term_mx = prhs[3];
        const mxArray *x_min_hard_term_mx = prhs[4];
        const mxArray *x_max_hard_term_mx = prhs[5];

        // declare
        ZonoCpp::ConZono Zc_x_term;
        Eigen::MatrixXd G_x_term, A_x_term, Pp_x_term;
        Eigen::VectorXd c_x_term, b_x_term, x_min_hard_term, x_max_hard_term;

        // get conzono
        getFullMatrixInput(mxGetField(Zc_x_term_mx, 0, "G"), G_x_term);
        getFullMatrixInput(mxGetField(Zc_x_term_mx, 0, "A"), A_x_term);
        getVectorInput(mxGetField(Zc_x_term_mx, 0, "c"), c_x_term);
        getVectorInput(mxGetField(Zc_x_term_mx, 0, "b"), b_x_term);
        Zc_x_term.set(G_x_term.sparseView(), c_x_term, A_x_term.sparseView(), b_x_term);

        // get Pp_x_term
        getFullMatrixInput(Pp_x_term_mx, Pp_x_term);

        // get x_min_hard_term and x_max_hard_term
        getVectorInput(x_min_hard_term_mx, x_min_hard_term);
        getVectorInput(x_max_hard_term_mx, x_max_hard_term);

        // call method
        instance->set_terminal_state_inequality_constraints(Zc_x_term, Pp_x_term.sparseView(), x_min_hard_term, x_max_hard_term);

        break;
    }
    case Action::Set_state_inequality_constraints_slack_vars_cost:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two input required: handle, state_ineq_slack_cost, sigma_max_x.");

        // get input data pointers
        const mxArray *state_ineq_slack_cost_mx = prhs[2];
        const mxArray *sigma_max_x_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> Qx_constraint_cost;
        std::vector<Eigen::Triplet<double>> Qx_constraint_cost_triplets;
        Qx_constraint_cost_triplets.reserve(getNNZ(state_ineq_slack_cost_mx));
        getSparseMatrixInput(state_ineq_slack_cost_mx, Qx_constraint_cost_triplets, &Qx_constraint_cost);

        Eigen::VectorXd sigma_max_x;
        getVectorInput(sigma_max_x_mx, sigma_max_x);

        // call method
        instance->set_state_inequality_constraints_slack_vars_cost(Qx_constraint_cost, sigma_max_x);

        break;
    }
    case Action::Set_input_inequality_constraints_slack_vars_cost:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two input required: handle, input_ineq_slack_cost, sigma_max_u.");

        // get input data pointers
        const mxArray *input_ineq_slack_cost_mx = prhs[2];
        const mxArray *sigma_max_u_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> Qu_constraint_cost;
        std::vector<Eigen::Triplet<double>> Qu_constraint_cost_triplets;
        Qu_constraint_cost_triplets.reserve(getNNZ(input_ineq_slack_cost_mx));
        getSparseMatrixInput(input_ineq_slack_cost_mx, Qu_constraint_cost_triplets, &Qu_constraint_cost);

        Eigen::VectorXd sigma_max_u;
        getVectorInput(sigma_max_u_mx, sigma_max_u);

        // call method
        instance->set_input_inequality_constraints_slack_vars_cost(Qu_constraint_cost, sigma_max_u);

        break;
    }
    case Action::Set_terminal_state_inequality_constraints_slack_vars_cost:
    {
        // check number of inputs
        if (nrhs != 2+2)
            mexErrMsgTxt("Two input required: handle, term_state_ineq_slack_cost, sigma_max_x_term.");

        // get input data pointers
        const mxArray *term_state_ineq_slack_cost_mx = prhs[2];
        const mxArray *sigma_max_x_term_mx = prhs[3];

        // get data
        Eigen::SparseMatrix<double> Qx_term_constraint_cost;
        std::vector<Eigen::Triplet<double>> Qx_term_constraint_cost_triplets;
        Qx_term_constraint_cost_triplets.reserve(getNNZ(term_state_ineq_slack_cost_mx));
        getSparseMatrixInput(term_state_ineq_slack_cost_mx, Qx_term_constraint_cost_triplets, &Qx_term_constraint_cost);

        Eigen::VectorXd sigma_max_x_term;
        getVectorInput(sigma_max_x_term_mx, sigma_max_x_term);

        // call method
        instance->set_terminal_state_inequality_constraints_slack_vars_cost(Qx_term_constraint_cost, sigma_max_x_term);

        break;
    }
    case Action::Set_QP_settings:
    {
        // check number of inputs
        if (nrhs != 2+1)
            mexErrMsgTxt("One input required: handle, QP_settings.");

        // get input data pointers
        const mxArray *qp_settings_ptr = prhs[2];

        // get data
        QP_IP_MPC::QP_settings qp_settings;
        qp_settings.iter_max = (unsigned int) mxGetScalar(mxGetField(qp_settings_ptr, 0, "iter_max"));
        qp_settings.mu_init = (double) mxGetScalar(mxGetField(qp_settings_ptr, 0, "mu_init"));
        qp_settings.mu_max = (double) mxGetScalar(mxGetField(qp_settings_ptr, 0, "mu_max"));
        qp_settings.mu_feas = (double) mxGetScalar(mxGetField(qp_settings_ptr, 0, "mu_feas"));
        qp_settings.mu_term = (double) mxGetScalar(mxGetField(qp_settings_ptr, 0, "mu_term"));
        qp_settings.gamma = (double) mxGetScalar(mxGetField(qp_settings_ptr, 0, "gamma"));
        qp_settings.t_ls = (double) mxGetScalar(mxGetField(qp_settings_ptr, 0, "t_ls"));
        qp_settings.T_max = (double) mxGetScalar(mxGetField(qp_settings_ptr, 0, "T_max"));

        // call method
        instance->set_qp_settings(qp_settings);

        break;
    }
    case Action::Set_MI_settings:
    {
        // check number of inputs
        if (nrhs != 2+1)
            mexErrMsgTxt("One input required: handle, MI_settings.");
        
        // get input data pointers
        const mxArray *mi_settings_ptr = prhs[2];

        // get data
        MI_QPIPMPC::MI_QPIPMPC_Settings mi_settings;
        mi_settings.max_iter_bb = (int) mxGetScalar(mxGetField(mi_settings_ptr, 0, "max_iter_bb"));
        mi_settings.conv_abs = (double) mxGetScalar(mxGetField(mi_settings_ptr, 0, "conv_abs"));
        mi_settings.conv_rel = (double) mxGetScalar(mxGetField(mi_settings_ptr, 0, "conv_rel"));
        mi_settings.eps_feas = (double) mxGetScalar(mxGetField(mi_settings_ptr, 0, "eps_feas"));
        mi_settings.verbose = (bool) mxGetScalar(mxGetField(mi_settings_ptr, 0, "verbose"));
        mi_settings.T_max = (double) mxGetScalar(mxGetField(mi_settings_ptr, 0, "T_max"));
        mi_settings.n_threads = (unsigned int) mxGetScalar(mxGetField(mi_settings_ptr, 0, "n_threads"));

        // set verbosity flag
        verbose = mi_settings.verbose;

        // call method
        instance->set_MI_settings(mi_settings);

        break;
    }
    case Action::Set_MPC_settings:
    {
        // check number of inputs
        if (nrhs != 2+1)
            mexErrMsgTxt("One input required: handle, MPC_settings.");

        // get input data pointers
        const mxArray *mpc_settings_ptr = prhs[2];

        // get data
        HybZonoMPC::MPC_Settings mpc_settings;
        mpc_settings.warm_start = (bool) mxGetScalar(mxGetField(mpc_settings_ptr, 0, "warm_start"));
        mpc_settings.u1_control = (bool) mxGetScalar(mxGetField(mpc_settings_ptr, 0, "u1_control"));
        mpc_settings.n_horizon = (int) mxGetScalar(mxGetField(mpc_settings_ptr, 0, "n_horizon"));

        // call method
        instance->set_mpc_settings(mpc_settings);

        break;
    }
    case Action::Set_distance_limit:
    {
        // check number of inputs
        if (nrhs != 2+1)
            mexErrMsgTxt("One input required: handle, dist_limit.");

        // get input data pointers
        const mxArray *dist_limit_mx = prhs[2];

        // get data
        double dist_limit = (double) mxGetScalar(dist_limit_mx);

        // call method
        instance->set_distance_limit(dist_limit);

        break;
    }
    case Action::Setup_hybzono_constraints:
    {
        mexErrMsgTxt("Not currently implemented.");
        break;
    }
    case Action::Setup_hybzono_Hrep_constraints:
    {
        // check number of inputs
        if (nrhs != 11+2)
            mexErrMsgTxt("11 inputs required for setup_hybzono_Hrep_constraints()");

        // declare inputs
        Eigen::SparseMatrix<double> Pp_3, Pc_3, Q_hz_3;
        Eigen::MatrixXd Pp_full_3, Pc_full_3;
        std::map<int, Eigen::MatrixXd> A_Hrep_3, A_Hrep_pos_3;
        std::map<int, Eigen::VectorXd> b_Hrep_3, b_Hrep_pos_3;
        std::vector<Eigen::Triplet<double>> Pp_triplets_3, Q_hz_triplets_3;
        ZonoCpp::HybZono X_hz_3;
        Eigen::MatrixXd Gc_3, Gb_3, Ac_3, Ab_3;
        Eigen::VectorXd c_3, b_3;
        double slack_min_hz_3=-1e4, slack_max_hz_3=1e4;
        Eigen::VectorXd cost_vec_3;

        // pointers
        const mxArray * X_hz_3_mx = prhs[2];
        const mxArray * A_Hrep_3_mx = prhs[3];
        const mxArray * b_Hrep_3_mx = prhs[4];
        const mxArray * Pc_3_mx = prhs[5];
        const mxArray * Q_hz_3_mx = prhs[6];
        const mxArray * slack_min_hz_3_mx = prhs[7];
        const mxArray * slack_max_hz_3_mx = prhs[8];
        const mxArray * cost_vec_3_mx = prhs[9];
        const mxArray * Pp_3_mx = prhs[10];
        const mxArray * A_Hrep_pos_3_mx = prhs[11];
        const mxArray * b_Hrep_pos_3_mx = prhs[12];

        // get hybzono
        getFullMatrixInput(mxGetField(X_hz_3_mx, 0, "Gc"), Gc_3);
        getFullMatrixInput(mxGetField(X_hz_3_mx, 0, "Gb"), Gb_3);
        getFullMatrixInput(mxGetField(X_hz_3_mx, 0, "Ac"), Ac_3);
        getFullMatrixInput(mxGetField(X_hz_3_mx, 0, "Ab"), Ab_3);
        getVectorInput(mxGetField(X_hz_3_mx, 0, "c"), c_3);
        getVectorInput(mxGetField(X_hz_3_mx, 0, "b"), b_3);
        X_hz_3.set(Gc_3.sparseView(), Gb_3.sparseView(), c_3, Ac_3.sparseView(), Ab_3.sparseView(), b_3); // {-1,1} binaries

        // get H-rep
        const mxArray *A_Hrep_mat_3_mxarray = mxGetField(A_Hrep_3_mx, 0, "A_Hrep_mat");
        const mxArray *A_Hrep_ind_3_mxarray = mxGetField(A_Hrep_3_mx, 0, "A_Hrep_ind");
        getFullMatrixMap(A_Hrep_mat_3_mxarray, A_Hrep_ind_3_mxarray, A_Hrep_3);

        const mxArray *b_Hrep_3_mat_mxarray = mxGetField(b_Hrep_3_mx, 0, "b_Hrep_mat");
        const mxArray *b_Hrep_3_ind_mxarray = mxGetField(b_Hrep_3_mx, 0, "b_Hrep_ind");
        getVectorMap(b_Hrep_3_mat_mxarray, b_Hrep_3_ind_mxarray, b_Hrep_3);

        // get Pc
        getFullMatrixInput(Pc_3_mx, Pc_full_3);
        Pc_3 = Pc_full_3.sparseView();

        // get slack cost
        if (!mxIsEmpty(Q_hz_3_mx))
        {
            Q_hz_triplets_3.reserve(getNNZ(Q_hz_3_mx));
            getSparseMatrixInput(Q_hz_3_mx, Q_hz_triplets_3, &Q_hz_3);
        }

        // get max slack
        if (!mxIsEmpty(slack_min_hz_3_mx))
            slack_min_hz_3 = (double) mxGetScalar(slack_min_hz_3_mx);

        if (!mxIsEmpty(slack_max_hz_3_mx))
            slack_max_hz_3 = (double) mxGetScalar(slack_max_hz_3_mx);

        // get cost vec
        if (!mxIsEmpty(cost_vec_3_mx))
            getVectorInput(cost_vec_3_mx, cost_vec_3);

        // get Pp
        if (!mxIsEmpty(Pp_3_mx))
        {
            getFullMatrixInput(Pp_3_mx, Pp_full_3);
            Pp_3 = Pp_full_3.sparseView();
        }

        // get H-rep pos
        const mxArray *A_Hrep_pos_mat_3_mxarray=nullptr, *A_Hrep_pos_ind_3_mxarray=nullptr;
        const mxArray *b_Hrep_pos_3_mat_mxarray=nullptr, *b_Hrep_pos_3_ind_mxarray=nullptr;
        if (!mxIsEmpty(A_Hrep_pos_3_mx) && !mxIsEmpty(b_Hrep_pos_3_mx))
        {
            A_Hrep_pos_mat_3_mxarray = mxGetField(A_Hrep_pos_3_mx, 0, "A_Hrep_pos_mat");
            A_Hrep_pos_ind_3_mxarray = mxGetField(A_Hrep_pos_3_mx, 0, "A_Hrep_pos_ind");
            getFullMatrixMap(A_Hrep_pos_mat_3_mxarray, A_Hrep_pos_ind_3_mxarray, A_Hrep_pos_3);

            b_Hrep_pos_3_mat_mxarray = mxGetField(b_Hrep_pos_3_mx, 0, "b_Hrep_pos_mat");
            b_Hrep_pos_3_ind_mxarray = mxGetField(b_Hrep_pos_3_mx, 0, "b_Hrep_pos_ind");
            getVectorMap(b_Hrep_pos_3_mat_mxarray, b_Hrep_pos_3_ind_mxarray, b_Hrep_pos_3);
        }

        // call method
        instance->setup_hybzono_Hrep_constraints(X_hz_3, A_Hrep_3, b_Hrep_3, Pc_3, Q_hz_3, slack_min_hz_3, slack_max_hz_3, cost_vec_3,
            Pp_3, A_Hrep_pos_3, b_Hrep_pos_3);

        break;
    }
    case Action::Set_region_dependent_disturbances:
    {
        // check number of inputs
        if (nrhs != 2+1)
            mexErrMsgTxt("One input required: handle, region_dependent_disturbances.");

        // get input data pointers
        const mxArray *region_dependent_disturbances_mx = prhs[2];

        // get data
        std::map<int, Eigen::VectorXd> d_regions;

        const mxArray *d_region_vec_mxarray = mxGetField(region_dependent_disturbances_mx, 0, "d_region_vec");
        const mxArray *d_region_ind_mxarray = mxGetField(region_dependent_disturbances_mx, 0, "d_region_ind");
        getVectorMap(d_region_vec_mxarray, d_region_ind_mxarray, d_regions);

        // call method
        instance->set_region_dependent_disturbances(d_regions);

        break;
    }
    case Action::Set_dist_between_regions:
    {
        // check number of inputs
        if (nrhs != 2+1)
            mexErrMsgTxt("One input required: handle, dist_between_regions.");

        // get input data pointers
        const mxArray *dist_between_regions_mx = prhs[2];

        // get input data pointersSet_terminal_hybzono_constraints
        // get data
        Eigen::MatrixXd d_mat;
        getFullMatrixInput(dist_between_regions_mx, d_mat);

        // call method
        instance->set_dist_between_regions(d_mat);

        break;
    }
    case Action::Set_terminal_hybzono_constraints:
    {
        // check number of inputs
        if (nrhs != 2+4)
            mexErrMsgTxt("Four inputs required: handle, X_hz_term, term_region_cost_vec, Pp_term, Q_term_slack_cost.");

        // get input data pointers
        const mxArray *X_hz_term_mx = prhs[2];
        const mxArray *term_region_cost_vec_mx = prhs[3];
        const mxArray *Pp_term_mx = prhs[4];
        const mxArray *Q_term_slack_cost_mx = prhs[5];

        // get hybzono
        ZonoCpp::HybZono X_hz_term;
        Eigen::MatrixXd Gc_term, Gb_term, Ac_term, Ab_term;
        Eigen::VectorXd c_term, b_term;
        getFullMatrixInput(mxGetField(X_hz_term_mx, 0, "Gc"), Gc_term);
        getFullMatrixInput(mxGetField(X_hz_term_mx, 0, "Gb"), Gb_term);
        getFullMatrixInput(mxGetField(X_hz_term_mx, 0, "Ac"), Ac_term);
        getFullMatrixInput(mxGetField(X_hz_term_mx, 0, "Ab"), Ab_term);
        getVectorInput(mxGetField(X_hz_term_mx, 0, "c"), c_term);
        getVectorInput(mxGetField(X_hz_term_mx, 0, "b"), b_term);
        X_hz_term.set(Gc_term.sparseView(), Gb_term.sparseView(), c_term, Ac_term.sparseView(), Ab_term.sparseView(), b_term); // {-1,1} binaries

        // get term_region_cost_vec
        Eigen::VectorXd term_region_cost_vec;
        getVectorInput(term_region_cost_vec_mx, term_region_cost_vec);

        // get Pp_term
        Eigen::SparseMatrix<double> Pp_term;
        Eigen::MatrixXd Pp_term_full;
        getFullMatrixInput(Pp_term_mx, Pp_term_full);
        Pp_term = Pp_term_full.sparseView();

        // get Q_term_slack_cost
        Eigen::SparseMatrix<double> Q_term_slack_cost;
        std::vector<Eigen::Triplet<double>> Q_term_slack_cost_triplets;
        Q_term_slack_cost_triplets.reserve(getNNZ(Q_term_slack_cost_mx));
        getSparseMatrixInput(Q_term_slack_cost_mx, Q_term_slack_cost_triplets, &Q_term_slack_cost);

        // call method
        instance->set_terminal_hybzono_constraints(X_hz_term, term_region_cost_vec, Pp_term, Q_term_slack_cost);

        break;
    }
    case Action::Build_controller:
    {
        try
        {
            // call method
            instance->build_controller();
        }
        catch (const std::exception &e)
        {
            mexErrMsgTxt(e.what());
        }

        break;
    }
    case Action::Control:
    {
        // check number of inputs
        if (nrhs != 2+2 && nrhs != 2+3)
            mexErrMsgTxt("Two or three inputs required: handle, x0, x_ref, u0");

        bool u0_specified = (nrhs == 2+3);

        // get input data pointers
        const mxArray *x0_mx = prhs[2];
        const mxArray *x_ref_mx = prhs[3];
        const mxArray *u0_mx = nullptr;
        if (u0_specified)
            u0_mx = prhs[4];

        // get x0
        Eigen::VectorXd x0;
        getVectorInput(x0_mx, x0);

        // get x_ref
        Eigen::MatrixXd x_ref_mat;
        getFullMatrixInput(x_ref_mx, x_ref_mat);

        std::vector<Eigen::VectorXd> x_ref;
        x_ref.reserve(x_ref_mat.cols());
        for (int i = 0; i < x_ref_mat.cols(); i++)
            x_ref.push_back(x_ref_mat.col(i));

        // u0
        Eigen::VectorXd u0;
        if (u0_specified)
            getVectorInput(u0_mx, u0);

        // call method
        std::pair<Eigen::VectorXd, bool> ctrl_out;
        if (u0_specified)
            ctrl_out = instance->control(x0, x_ref, u0);
        else
            ctrl_out = instance->control(x0, x_ref);
        
        Eigen::VectorXd u_out = ctrl_out.first;
        bool is_feasible_out = ctrl_out.second;

        // create output
        plhs[0] = mxCreateDoubleMatrix(u_out.size(), 1, mxREAL);
        copyEigenVector2DblPtr(u_out, mxGetPr(plhs[0]));

        // optional flags
        plhs[1] = mxCreateDoubleScalar((double) is_feasible_out);

        // verbose output if applicable
        if (verbose)
        {
            verbose_output = instance->get_verbose_output();
            mexPrintf("%s\n", verbose_output.c_str());
        }

        break;
    }
    case Action::Control_LTV:
    {
        // check number of inputs
        if (nrhs != 5+2)
            mexErrMsgTxt("Five inputs required: handle, x0, x_ref, u0, A_dyn_vec, B_dyn_vec.");

        // get input data pointers
        const mxArray *x0_ltv_mx = prhs[2];
        const mxArray *x_ref_ltv_mx = prhs[3];
        const mxArray *u0_ltv_mx = prhs[4];
        const mxArray *A_dyn_vec_control_mx = prhs[5];
        const mxArray *B_dyn_vec_control_mx = prhs[6];

        // get x0
        Eigen::VectorXd x0_ltv;
        getVectorInput(x0_ltv_mx, x0_ltv);

        // get x_ref
        Eigen::MatrixXd x_ref_ltv_mat;
        getFullMatrixInput(x_ref_ltv_mx, x_ref_ltv_mat);

        std::vector<Eigen::VectorXd> x_ref_ltv;
        x_ref_ltv.reserve(x_ref_ltv_mat.cols());
        for (int i = 0; i < x_ref_ltv_mat.cols(); i++)
            x_ref_ltv.push_back(x_ref_ltv_mat.col(i));

        // get x0
        Eigen::VectorXd u0_ltv;
        getVectorInput(u0_ltv_mx, u0_ltv);

        // get A_dyn, B_dyn
        std::vector<Eigen::SparseMatrix<double>> A_dyn_vec_control, B_dyn_vec_control;
        std::vector<Eigen::Triplet<double>> A_dyn_vec_control_triplets, B_dyn_vec_control_triplets;
        Eigen::SparseMatrix<double> A_dyn_control_el, B_dyn_control_el;

        // loop through cell array to get matrices
        for (int i = 0; i < mxGetNumberOfElements(A_dyn_vec_control_mx); i++)
        {
            A_dyn_vec_control_triplets.clear();
            B_dyn_vec_control_triplets.clear();
            const mxArray *A_dyn_control_mx = mxGetCell(A_dyn_vec_control_mx, i);
            const mxArray *B_dyn_control_mx = mxGetCell(B_dyn_vec_control_mx, i);
            getSparseMatrixInput(A_dyn_control_mx, A_dyn_vec_control_triplets, &A_dyn_control_el);
            getSparseMatrixInput(B_dyn_control_mx, B_dyn_vec_control_triplets, &B_dyn_control_el);
            A_dyn_vec_control.push_back(A_dyn_control_el);
            B_dyn_vec_control.push_back(B_dyn_control_el);
        }

        // call method
        std::pair<Eigen::VectorXd, bool> ctrl_ltv_out = instance->control(x0_ltv, x_ref_ltv, u0_ltv, A_dyn_vec_control, B_dyn_vec_control);
        Eigen::VectorXd u_ltv_out = ctrl_ltv_out.first;
        bool is_feasible_ltv_out = ctrl_ltv_out.second;

        // create output
        plhs[0] = mxCreateDoubleMatrix(u_ltv_out.size(), 1, mxREAL);
        copyEigenVector2DblPtr(u_ltv_out, mxGetPr(plhs[0]));

        // optional flags
        plhs[1] = mxCreateDoubleScalar((double) is_feasible_ltv_out);

        if (verbose)
        {
            verbose_output = instance->get_verbose_output();
            mexPrintf("%s\n", verbose_output.c_str());
        }

        break;
    }
    case Action::Get_trajectory:
    {
        // call method
        std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>> traj = instance->get_trajectory();

        // convert into dense matrices
        Eigen::MatrixXd x_traj(traj.first[0].size(), traj.first.size());
        for (int i=0; i<traj.first.size(); i++)
            x_traj.col(i) = traj.first[i];

        Eigen::MatrixXd u_traj(traj.second[0].size(), traj.second.size());
        for (int i=0; i<traj.second.size(); i++)
            u_traj.col(i) = traj.second[i];

        // create output
        plhs[0] = mxCreateDoubleMatrix(x_traj.rows(), x_traj.cols(), mxREAL);
        plhs[1] = mxCreateDoubleMatrix(u_traj.rows(), u_traj.cols(), mxREAL);
        copyEigenMatrix2DblPtr(x_traj, mxGetPr(plhs[0]));
        copyEigenMatrix2DblPtr(u_traj, mxGetPr(plhs[1]));

        break;
    }
    case Action::Get_terminal_region_selection:
    {
        // call method
        int term_region_selection = instance->get_terminal_region_selection();

        // return value
        plhs[0] = mxCreateDoubleScalar((double) term_region_selection);

        break;
    }
    case Action::Get_miqp_results:
    {
        // call method
        MI_QPIPMPC::Results miqp_results = instance->get_miqp_results();

        // create output struct
        int size_results_fields = sizeof(MI_RESULTS_FIELDS)/sizeof(MI_RESULTS_FIELDS[0]);
        plhs[0] = mxCreateStructMatrix(1, 1, size_results_fields, MI_RESULTS_FIELDS);

        mxSetField(plhs[0], 0, "upper_glob", mxCreateDoubleScalar(miqp_results.upper_glob));
        mxSetField(plhs[0], 0, "lower_glob", mxCreateDoubleScalar(miqp_results.lower_glob));
        mxSetField(plhs[0], 0, "run_time", mxCreateDoubleScalar(miqp_results.run_time));
        mxSetField(plhs[0], 0, "status", mxCreateDoubleScalar(mi_status_2_double(miqp_results.status)));
        mxSetField(plhs[0], 0, "qp_solve_time", mxCreateDoubleScalar(miqp_results.qp_solve_time));
        mxSetField(plhs[0], 0, "qp_iter_avg", mxCreateDoubleScalar((double) miqp_results.qp_iter_avg));
        mxSetField(plhs[0], 0, "iter_num", mxCreateDoubleScalar((double) miqp_results.iter_num));

        // set x
        plhs[1] = mxCreateDoubleMatrix(miqp_results.x.rows(), 1, mxREAL);
        copyEigenVector2DblPtr(miqp_results.x, mxGetPr(plhs[1]));

        // set region_vec
        plhs[2] = mxCreateDoubleMatrix(miqp_results.region_vec.size(), 1, mxREAL);
        copyStdVector2DblPtr(miqp_results.region_vec, mxGetPr(plhs[2]));

        break;
    }
    default:
    {
        mexErrMsgTxt(("Unhandled action: " + actionStr).c_str());
        break;
    }
}

}

handle_type getHandle(int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2 || mxGetNumberOfElements(prhs[1]) != 1) // mxIsScalar in R2015a+
        mexErrMsgTxt("Specify an instance with an integer handle.");
    return static_cast<handle_type>(mxGetScalar(prhs[1]));
}

instanceMap_type::const_iterator checkHandle(const instanceMap_type& m, handle_type h)
{
    auto it = m.find(h);

    if (it == m.end()) {
        std::stringstream ss; ss << "No instance corresponding to handle " << h << " found.";
        mexErrMsgTxt(ss.str().c_str());
    }

    return it;
}
