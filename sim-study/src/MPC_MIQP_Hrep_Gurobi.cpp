#include "MPC_MIQP_Hrep_Gurobi.hpp"

using namespace HybZonoMPC;

// destructor
MPC_MIQP_Hrep_Gurobi::~MPC_MIQP_Hrep_Gurobi()
{
    delete [] x_gur;
    delete model_ptr;
}

// setup
void MPC_MIQP_Hrep_Gurobi::set_MI_settings(const Gurobi_Settings &settings)
{
    gurobi_settings = settings;
}

void MPC_MIQP_Hrep_Gurobi::build_controller()
{
    MPC_MIQP_Hrep::build_controller();
    make_full_system_matrices();
    set_gurobi_settings();
}

void MPC_MIQP_Hrep_Gurobi::set_gurobi_settings()
{
    // set number of threads
    env.set(GRB_IntParam_Threads, gurobi_settings.n_threads);

    // set time limit
    if (gurobi_settings.T_max == 0)
        env.set(GRB_DoubleParam_TimeLimit, GRB_INFINITY);
    else
        env.set(GRB_DoubleParam_TimeLimit, gurobi_settings.T_max);

    // set verbosity
    if (gurobi_settings.verbose)
        env.set(GRB_IntParam_OutputFlag, 1);
    else
        env.set(GRB_IntParam_OutputFlag, 0);

    // set MIP gap
    env.set(GRB_DoubleParam_MIPGap, gurobi_settings.conv_rel);
    env.set(GRB_DoubleParam_MIPGapAbs, gurobi_settings.conv_abs);
}

void MPC_MIQP_Hrep_Gurobi::set_terminal_hybzono_constraints(const ZonoCpp::HybZono &X_hz_term, const Eigen::Ref<const Eigen::VectorXd> term_region_cost_vec,
            const Eigen::SparseMatrix<double> &Pp_term, const Eigen::SparseMatrix<double> &Q_term_slack_cost)
{
    throw std::runtime_error("Terminal hybzono constraints not currently supported using Gurobi as MIQP solver");
}

int MPC_MIQP_Hrep_Gurobi::get_terminal_region_selection()
{
    throw std::runtime_error("Terminal hybzono constraints not currently supported using Gurobi as MIQP solver");
}

// control methods (LTI)
std::pair<Eigen::VectorXd, bool> MPC_MIQP_Hrep_Gurobi::control(const Eigen::Ref<const Eigen::VectorXd> x, const std::vector<Eigen::VectorXd> &x_ref)
{
    // set x0 and x_ref
    this->x0 = x;
    this->x_ref = x_ref;

    // constraint update
    make_IC_constraint();

    // gradient update
    make_cost_gradient_from_reference();

    // LTV update
    if (LTV)
    {
        update_equality_constraints_LTV();
    }

    // solve
    make_optim_problem();
    bool feasible = solve_optimization_problem();

    // return control input for time zero (or 1 if using u1 control)
    if (u1_control)
        return std::make_pair(u_vec[1], feasible);
    else
        return std::make_pair(u_vec[0], feasible);
}

std::pair<Eigen::VectorXd, bool>  MPC_MIQP_Hrep_Gurobi::control(const Eigen::Ref<const Eigen::VectorXd> x)
{
    // set x_ref
    this->x_ref = std::vector<Eigen::VectorXd>(n_horizon, Eigen::VectorXd::Zero(n_state));

    // call control method
    return control(x, this->x_ref);
}

// control methods (LTV)
std::pair<Eigen::VectorXd, bool>  MPC_MIQP_Hrep_Gurobi::control(const Eigen::Ref<const Eigen::VectorXd> x, const std::vector<Eigen::VectorXd> &x_ref, 
            const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec, const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec)
{
    // set LTV flag
    LTV = true;

    // updated dynamics matrices
    this->A_dyn_vec = A_dyn_vec;
    this->B_dyn_vec = B_dyn_vec;
    this->A_dyn = *(A_dyn_vec.begin());
    this->B_dyn = *(B_dyn_vec.begin());

    // call control method
    return control(x, x_ref);
}

std::pair<Eigen::VectorXd, bool>  MPC_MIQP_Hrep_Gurobi::control(const Eigen::Ref<const Eigen::VectorXd> x,
            const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec, const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec)
{
    // set LTV flag
    LTV = true;

    // updated dynamics matrices
    this->A_dyn_vec = A_dyn_vec;
    this->B_dyn_vec = B_dyn_vec;
    this->A_dyn = *(A_dyn_vec.begin());
    this->B_dyn = *(B_dyn_vec.begin());

    // set x_ref
    this->x_ref = std::vector<Eigen::VectorXd>(n_horizon, Eigen::VectorXd::Zero(n_state));

    // call control method
    return control(x, this->x_ref);
}

void MPC_MIQP_Hrep_Gurobi::make_full_system_matrices()
{
    // declare
    std::vector<Eigen::Triplet<double>> tripvec;
    int i_offset, j_offset;

    // P
    Eigen::SparseMatrix<double> * P_ptr = nullptr;
    tripvec.clear();
    i_offset = 0;
    j_offset = 0;
    for (int k=0; k<=n_horizon; k++)
    {   
        if (P_i_vec.count(k))
            P_ptr = &P_i_vec.at(k);
        else
            P_ptr = &P_i_nom;
        get_triplets_offset(*P_ptr, tripvec, i_offset, j_offset);
        i_offset += P_ptr->rows();
        j_offset += P_ptr->cols();
    }
    P.resize(n_y_tot, n_y_tot);
    P.setFromTriplets(tripvec.begin(), tripvec.end());

    // q
    Eigen::VectorXd * q_ptr = nullptr;
    q.resize(n_y_tot);
    i_offset = 0;
    for (int k=0; k<=n_horizon; k++)
    {
        if (q_i_vec.count(k))
            q_ptr = &q_i_vec.at(k);
        else
            q_ptr = &q_i_nom;
        q.segment(i_offset, q_ptr->size()) = *q_ptr;
        i_offset += q_ptr->size();
    }

    // A
    Eigen::SparseMatrix<double> * C_ptr = nullptr;
    Eigen::SparseMatrix<double> * D_ptr = nullptr;
    tripvec.clear();
    i_offset = 0;
    j_offset = 0;
    for (int k=0; k<n_horizon; k++)
    {
        if (C_i_vec.count(k))
            C_ptr = &C_i_vec.at(k);
        else
            C_ptr = &C_i_nom;
        if (D_i_vec.count(k+1))
            D_ptr = &D_i_vec.at(k+1);
        else
            D_ptr = &D_i_nom;

        get_triplets_offset(*C_ptr, tripvec, i_offset, j_offset);
        get_triplets_offset(*D_ptr, tripvec, i_offset, j_offset+C_ptr->cols());

        i_offset += C_ptr->rows();
        j_offset += C_ptr->cols();
    }
    A.resize(n_eq_tot, n_y_tot);
    A.setFromTriplets(tripvec.begin(), tripvec.end());

    // b
    Eigen::VectorXd * crhs_ptr = nullptr;
    b.resize(n_eq_tot);
    i_offset = 0;
    for (int k=0; k<n_horizon; k++)
    {
        if (crhs_i_vec.count(k))
            crhs_ptr = &crhs_i_vec.at(k);
        else
            crhs_ptr = &crhs_i_nom;
        b.segment(i_offset, crhs_ptr->size()) = -1*(*crhs_ptr);
        i_offset += crhs_ptr->size();
    }

    // G
    Eigen::SparseMatrix<double> * G_ptr = nullptr;
    tripvec.clear();
    i_offset = 0;
    j_offset = 0;
    for (int k=0; k<=n_horizon; k++)
    {
        if (G_i_vec.count(k))
            G_ptr = &G_i_vec.at(k);
        else
            G_ptr = &G_i_nom;
        get_triplets_offset(*G_ptr, tripvec, i_offset, j_offset);
        i_offset += G_ptr->rows();
        j_offset += G_ptr->cols();
    }
    G.resize(n_ineq_tot, n_y_tot);
    G.setFromTriplets(tripvec.begin(), tripvec.end());

    // w
    Eigen::VectorXd * w_ptr = nullptr;
    w.resize(n_ineq_tot);
    i_offset = 0;
    for (int k=0; k<=n_horizon; k++)
    {
        if (w_i_vec.count(k))
            w_ptr = &w_i_vec.at(k);
        else
            w_ptr = &w_i_nom;
        w.segment(i_offset, w_ptr->size()) = *w_ptr;
        i_offset += w_ptr->size();
    }
}

void MPC_MIQP_Hrep_Gurobi::make_optim_problem()
{
    // compute system matrices
    make_full_system_matrices();

    // create model
    delete model_ptr;
    model_ptr = new GRBModel(env);

    // create optimization variables
    double * x_lb = new double [n_y_tot];
    double * x_ub = new double [n_y_tot];
    char * x_type = new char [n_y_tot];
    for (int i=0; i<n_y_tot; i++)
    {
        x_lb[i] = -1*GRB_INFINITY;
        x_ub[i] = GRB_INFINITY;
        x_type[i] = GRB_CONTINUOUS;
    }

    // declare certain variables as binary
    for (int k=0; k<=n_horizon; k++)
    {
        for (int i : idx_binvar.at(k))
        {
            x_type[i] = GRB_BINARY;
            x_lb[i] = 0;
            x_ub[i] = 1;
        }
    }

    // Gurobi variables
    delete [] x_gur;
    x_gur = model_ptr->addVars(x_lb, x_ub, NULL, x_type, NULL, n_y_tot);

    // objective
    quad_cost(model_ptr, x_gur, &P, q);

    // constraints
    lin_cons(model_ptr, x_gur, Eigen::MatrixXd(A), b, Eigen::MatrixXd(G), w);

    // free memory
    delete [] x_lb;
    delete [] x_ub;
    delete [] x_type;
}

bool MPC_MIQP_Hrep_Gurobi::solve_optimization_problem()
{
    // solve and time
    auto start = std::chrono::high_resolution_clock::now();
    model_ptr->optimize();
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    miqp_results.run_time = 1e-6 * ((double) elapsed.count());
    miqp_results.iter_num = (int) model_ptr->get(GRB_DoubleAttr_NodeCount);

    // status handling
    std::stringstream ss; // declare
    switch (model_ptr->get(GRB_IntAttr_Status))
    {
        case GRB_OPTIMAL:
        {
            miqp_results.status = MI_QPIPMPC::MIQP_Status::MIQP_SOLVED;
            break;
        }
        case GRB_INFEASIBLE: case GRB_INF_OR_UNBD: case GRB_UNBOUNDED:
        {
            miqp_results.status = MI_QPIPMPC::MIQP_Status::MIQP_INFEASIBLE;
            break;
        }
        case GRB_LOADED:
        {
            miqp_results.status = MI_QPIPMPC::MIQP_Status::MIQP_NO_SOL;
            throw std::runtime_error("Gurobi optimization not solved");
            break;
        }
        case GRB_CUTOFF: case GRB_ITERATION_LIMIT: case GRB_TIME_LIMIT: case GRB_NODE_LIMIT: case GRB_SOLUTION_LIMIT:
            case GRB_NUMERIC: case GRB_SUBOPTIMAL: case GRB_USER_OBJ_LIMIT: case GRB_WORK_LIMIT: case GRB_MEM_LIMIT:
        {
            miqp_results.status = MI_QPIPMPC::MIQP_Status::MIQP_MAX_ITER_FEASIBLE;
            break;
        }
        case GRB_INTERRUPTED: case GRB_INPROGRESS:
        {
            throw std::runtime_error("Gurobi optimization interrupted");
            break;
        }
        default:
        {
            ss << "Unknown Gurobi optimization status: " << model_ptr->get(GRB_IntAttr_Status);
            throw std::runtime_error(ss.str());
            break;
        }
    }

    // check if feasible
    bool feasible = miqp_results.status != MI_QPIPMPC::MIQP_Status::MIQP_INFEASIBLE;

    // get solution vector
    solution.resize(n_y_tot);
    if (feasible)
    {
        for (int i=0; i<n_y_tot; i++)
            solution(i) = x_gur[i].get(GRB_DoubleAttr_X);
    }
    else
    {
        solution = Eigen::VectorXd::Zero(n_y_tot);
    }
    miqp_results.x = solution;

    // get region vector
    miqp_results.region_vec.clear();
    for (int k=0; k<=n_horizon; k++)
    {
        Eigen::VectorXd regions_k = solution(idx_binvar.at(k));
        double max_val = regions_k.maxCoeff();
        int region = -1;
        for (int i=0; i<regions_k.size(); i++)
        {
            if (regions_k(i) == max_val)
            {
                region = i+1; // one-based indexing
                break;
            }
        }
        miqp_results.region_vec.push_back(region);
    }

    // get objective
    if (feasible)
        miqp_results.upper_glob = model_ptr->get(GRB_DoubleAttr_ObjVal);
    else
        miqp_results.upper_glob = QP::inf;

    // infeasibility handling
    get_state_input_trajectories();

    // unused results struct parameters
    miqp_results.qp_iter_avg = 0;
    miqp_results.qp_solve_time = 0;

    // update u0 if using u1 control
    if (u1_control)
        u0 = u_vec[1];

    // return
    return feasible;
}

// utilities
void MPC_MIQP_Hrep_Gurobi::quad_cost(GRBModel * model_ptr, const GRBVar *x, 
    const Eigen::SparseMatrix<double> * P, 
    const Eigen::Ref<const Eigen::VectorXd> q) 
{
    // init cost
    GRBQuadExpr obj = 0;

    // quadratic cost
    for (int i = 0; i < P->outerSize(); i++) 
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(*P, i); it; ++it) 
        {
            obj += 0.5 * it.value() * x[it.row()] * x[it.col()];
        }
    }

    // linear cost
    for (int i = 0; i < q.size(); i++) 
    {
        obj += q(i) * x[i];
    }

    // set objective
    model_ptr->setObjective(obj);
}

void MPC_MIQP_Hrep_Gurobi::lin_cons(GRBModel * model_ptr, const GRBVar *x, 
    const Eigen::Ref<const Eigen::MatrixXd> A_eq, const Eigen::Ref<const Eigen::VectorXd> b_eq,
    const Eigen::Ref<const Eigen::MatrixXd> A_ineq, const Eigen::Ref<const Eigen::VectorXd> b_ineq) 
{
    // init constraint expression
    GRBLinExpr expr;

    // add inequality constraints
    for (int i=0; i<A_ineq.rows(); i++) 
    {
        expr = 0;
        for (int j=0; j<A_ineq.cols(); j++) {
            expr += A_ineq(i,j) * x[j];
        }
        model_ptr->addConstr(expr <= b_ineq(i));
    }

    // add equality constraints
    for (int i=0; i<A_eq.rows(); i++) 
    {
        expr = 0;
        for (int j=0; j<A_eq.cols(); j++) {
            expr += A_eq(i,j) * x[j];
        }
        model_ptr->addConstr(expr == b_eq(i));
    }
}