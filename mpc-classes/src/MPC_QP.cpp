#include "MPC_QP.hpp"

using namespace HybZonoMPC;

// set methods
void MPC_QP::set_dyn_matrices(const Eigen::SparseMatrix<double> &A_dyn, const Eigen::SparseMatrix<double> &B_dyn)
{
    // set dynamics matrices
    this->A_dyn = A_dyn;
    this->B_dyn = B_dyn;

    // problem dimensions
    this->n_state = this->A_dyn.rows();
    this->n_input = this->B_dyn.cols();
}

void MPC_QP::set_dyn_matrices(const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec, const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec)
{
    // set flag
    this->LTV = true;

    // set dynamics matrices
    this->A_dyn_vec = A_dyn_vec;
    this->B_dyn_vec = B_dyn_vec;
    this->A_dyn = *(A_dyn_vec.begin());
    this->B_dyn = *(B_dyn_vec.begin());
    
    // problem dimensions
    this->n_state = this->A_dyn.rows();
    this->n_input = this->B_dyn.cols();
}

void MPC_QP::set_stage_cost(const Eigen::SparseMatrix<double> &Q_cost, const Eigen::SparseMatrix<double> &R_cost, 
        const Eigen::SparseMatrix<double> &N_cost, const Eigen::Ref<const Eigen::VectorXd> q_x_cost,
        const Eigen::Ref<const Eigen::VectorXd> q_u_cost)
{
    // set stage cost matrices and vectors
    // J_k = x'*Q*x + u'*R*u + 2*x'*N*u + 2*q_x'*x + 2*q_u'*u
    this->Q_cost = Q_cost;
    this->R_cost = R_cost;

    if (N_cost.rows() == 0)
        this->N_cost = Eigen::SparseMatrix<double>(this->Q_cost.rows(), this->R_cost.rows());
    else
        this->N_cost = N_cost;

    if (q_x_cost.size() == 0)
        this->q_x_cost = Eigen::VectorXd::Zero(this->Q_cost.rows());
    else
        this->q_x_cost = q_x_cost;

    if (q_u_cost.size() == 0)
        this->q_u_cost = Eigen::VectorXd::Zero(this->R_cost.rows());
    else
        this->q_u_cost = Eigen::VectorXd::Zero(this->R_cost.rows());

    // if terminal cost not already specified, use Q matrix
    if (!terminal_cost_specified)
    {
        this->P_cost = Q_cost;
        this->q_xN_cost = q_x_cost;
    }
}

void MPC_QP::set_terminal_cost(const Eigen::SparseMatrix<double> &P_cost, const Eigen::Ref<const Eigen::VectorXd> q_xN_cost)
{
    // set terminal cost matrices and vectors
    // J_N = x_N'*P*x_N + 2*q_x_N'*x
    this->P_cost = P_cost;

    if (q_xN_cost.size() == 0)
        this->q_xN_cost = Eigen::VectorXd::Zero(this->P_cost.rows());
    else
        this->q_xN_cost = q_xN_cost;

    this->terminal_cost_specified = true;
}

void MPC_QP::set_state_inequality_constraints(const Eigen::SparseMatrix<double> &Ax_ineq, const Eigen::VectorXd &bx_ineq)
{
    // set state inequality constraints
    this->Ax_ineq = Ax_ineq;
    this->bx_ineq = bx_ineq;

    // if terminal constraint not already specified, use these matrices
    if (!terminal_state_ineq_specified)
    {
        this->Ax_term_ineq = Ax_ineq;
        this->bx_term_ineq = bx_ineq;
    }

    // set flag
    this->state_ineq_specified = true;
    this->Hrep_constraints = true;
}

void MPC_QP::set_input_inequality_constraints(const Eigen::SparseMatrix<double> &Au_ineq, const Eigen::VectorXd &bu_ineq)
{
    // set input inequality constraints
    this->Au_ineq = Au_ineq;
    this->bu_ineq = bu_ineq;

    // set flag
    this->input_ineq_specified = true;
    this->Hrep_constraints = true;
}

void MPC_QP::set_terminal_state_inequality_constraints(const Eigen::SparseMatrix<double> &Ax_term_ineq, const Eigen::Ref<const Eigen::VectorXd> bx_term_ineq)
{
    // set terminal state inequality constraints
    this->Ax_term_ineq = Ax_term_ineq;
    this->bx_term_ineq = bx_term_ineq;

    // set flag
    this->terminal_state_ineq_specified = true;
    this->Hrep_constraints = true;
}

void MPC_QP::set_state_inequality_constraints(const ZonoCpp::ConZono &Zc_x, const Eigen::SparseMatrix<double> &Pp_x,
    const Eigen::Ref<const Eigen::VectorXd> x_min_hard, const Eigen::Ref<const Eigen::VectorXd> x_max_hard)
{
    // set state inequality constraints
    this->Zc_x = Zc_x;
    if (!this->Zc_x.zero_one_generators)
        this->Zc_x.convert_generator_range();

    this->Pp_x = Pp_x;
    this->x_min_hard = x_min_hard;
    this->x_max_hard = x_max_hard;

    // set flag
    this->state_ineq_specified = true;
    this->conzono_constraints = true;
}

void MPC_QP::set_input_inequality_constraints(const ZonoCpp::ConZono &Zc_u, const Eigen::SparseMatrix<double> &Pp_u,
    const Eigen::Ref<const Eigen::VectorXd> u_min_hard, const Eigen::Ref<const Eigen::VectorXd> u_max_hard)
{
    // set input inequality constraints
    this->Zc_u = Zc_u;
    if (!this->Zc_u.zero_one_generators)
        this->Zc_u.convert_generator_range();

    this->Pp_u = Pp_u;
    this->u_min_hard = u_min_hard;
    this->u_max_hard = u_max_hard;

    // set flag
    this->input_ineq_specified = true;
    this->conzono_constraints = true;
}

void MPC_QP::set_terminal_state_inequality_constraints(const ZonoCpp::ConZono &Zc_x_term, const Eigen::SparseMatrix<double> &Pp_x_term,
    const Eigen::Ref<const Eigen::VectorXd> x_term_min_hard, const Eigen::Ref<const Eigen::VectorXd> x_term_max_hard)
{
    // set terminal state inequality constraints
    this->Zc_x_term = Zc_x_term;
    if (!this->Zc_x_term.zero_one_generators)
        this->Zc_x_term.convert_generator_range();

    this->Pp_x_term = Pp_x_term;
    this->x_term_min_hard = x_term_min_hard;
    this->x_term_max_hard = x_term_max_hard;

    // set flag
    this->terminal_state_ineq_specified = true;
    this->conzono_constraints = true;
}

void MPC_QP::set_state_inequality_constraints_slack_vars_cost(const Eigen::SparseMatrix<double> &Qx_ineq_cost, Eigen::Ref<const Eigen::VectorXd> sigma_max_x)
{
    // set state inequality constraints slack variables cost
    this->Qx_ineq_cost = Qx_ineq_cost;
    this->sigma_max_x = sigma_max_x;

    // if terminal constraint softening cost is not specified, use this
    if (!soft_term_state_ineq)
    {
        this->Qx_term_ineq_cost = Qx_ineq_cost;
        this->sigma_max_x_term = sigma_max_x;
        this->soft_term_state_ineq = true;
    }

    // set flag
    this->soft_state_ineq = true;
}

void MPC_QP::set_input_inequality_constraints_slack_vars_cost(const Eigen::SparseMatrix<double> &Qu_ineq_cost, Eigen::Ref<const Eigen::VectorXd> sigma_max_u)
{
    // set input inequality constraints slack variables cost
    this->Qu_ineq_cost = Qu_ineq_cost;
    this->sigma_max_u = sigma_max_u;

    // set flag
    this->soft_input_ineq = true;
}

void MPC_QP::set_terminal_state_inequality_constraints_slack_vars_cost(const Eigen::SparseMatrix<double> &Qx_term_ineq_cost, Eigen::Ref<const Eigen::VectorXd> sigma_max_x_term)
{
    // set terminal state inequality constraints slack variables cost
    this->Qx_term_ineq_cost = Qx_term_ineq_cost;
    this->sigma_max_x_term = sigma_max_x_term;

    // set flag
    this->soft_term_state_ineq = true;
}

void MPC_QP::set_mpc_settings(const MPC_Settings &mpc_settings)
{
    // set mpc settings
    this->n_horizon = mpc_settings.n_horizon;
    this->u1_control = mpc_settings.u1_control;
    this->ref_dep_term_constraint = mpc_settings.ref_dep_term_constraint;
}

void MPC_QP::set_qp_settings(const QP_IP_MPC::QP_settings &qp_settings)
{
    // set qp settings
    this->qp_settings = qp_settings;
}

// build controller
void MPC_QP::build_controller()
{
    // check for constraint compatability
    if (conzono_constraints && Hrep_constraints)
    {
        throw std::invalid_argument("Mixing H-rep and conzono constraints is not supported");
        return;
    }
    if (!conzono_constraints && !Hrep_constraints)
    {
        conzono_constraints = true;
    }
    if (ref_dep_term_constraint)
    {
        if (!conzono_constraints)
        {
            throw std::invalid_argument("Reference dependent terminal constraints require conzono constraints");
            return;
        }
        ZonoCpp::ConZono * Zc_x_term_ptr = terminal_state_ineq_specified ? &Zc_x_term : &Zc_x;
        if (Zc_x_term_ptr->n != this->n_state)
        {
            throw std::invalid_argument("Terminal state inequality constraints must have same dimension as state when using reference dependent terminal constraints");
            return;
        }
    }

    // if constraints aren't specified, build them up as +/- inf
    if (Hrep_constraints)
    {
        if (!state_ineq_specified)
        {
            Ax_ineq = Eigen::SparseMatrix<double>(1, n_state);
            bx_ineq = QP::inf*Eigen::VectorXd::Ones(1);
        }
        if (!terminal_state_ineq_specified)
        {
            Ax_term_ineq = Eigen::SparseMatrix<double>(1, n_state);
            bx_term_ineq = QP::inf*Eigen::VectorXd::Ones(1);
        }
        if (!input_ineq_specified)
        {
            Au_ineq = Eigen::SparseMatrix<double>(1, n_input);
            bu_ineq = QP::inf*Eigen::VectorXd::Ones(1);
        }
    }
    
    // initialize x0 and x_ref
    if (x0.size() != n_state)
        x0 = Eigen::VectorXd::Zero(n_state);
    if (x_ref.size() != n_horizon)
        x_ref = std::vector<Eigen::VectorXd>(n_horizon, Eigen::VectorXd::Zero(n_state));
    if (u0.size() != n_input)
        u0 = Eigen::VectorXd::Zero(n_input);

    // check problem validity
    auto [valid, msg] = check_problem_validity();
    if (!valid)
    {
        throw std::invalid_argument("Problem not valid: " + msg);
        return;
    }

    // create matrices for optimization problem
    make_inequality_constraints();
    make_equality_constraints();
    make_cost_function();

    // get state dimensions at each time step
    get_problem_dimensions();

    // initial solver setup and configuration
    configure_solver();
}

// control methods (LTI)
std::pair<Eigen::VectorXd, bool> MPC_QP::control(const Eigen::Ref<const Eigen::VectorXd> x, 
    const std::vector<Eigen::VectorXd> &x_ref, const Eigen::Ref<const Eigen::VectorXd> u,
    const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec, 
    const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec)
{
    // set x0
    this->x0 = x;

    // handle optional inputs
    if (x_ref.size() == 0)
        this->x_ref = std::vector<Eigen::VectorXd>(n_horizon, Eigen::VectorXd::Zero(n_state));
    else if (x_ref.size() == 1)
        this->x_ref = std::vector<Eigen::VectorXd>(n_horizon, x_ref[0]);
    else if (x_ref.size() == n_horizon)
        this->x_ref = x_ref;
    else
        throw std::invalid_argument("Invalid reference trajectory length");

    if (u.size() == n_input && u1_control)
        this->u0 = u;
    else if (u1_control)
        throw std::invalid_argument("Require valid control input be specified");
    else if (u.size() != 0)
        throw std::invalid_argument("Invalid control input size");

    if (A_dyn_vec.size() == n_horizon && B_dyn_vec.size() == n_horizon)
    {
        // set LTV flag
        this->LTV = true;

        // updated dynamics matrices
        this->A_dyn_vec = A_dyn_vec;
        this->B_dyn_vec = B_dyn_vec;
        this->A_dyn = *(A_dyn_vec.begin());
        this->B_dyn = *(B_dyn_vec.begin());   
    }

    // constraint updates
    make_IC_constraint();
    if (ref_dep_term_constraint)
        update_ref_dep_term_constraint();

    // gradient update
    make_cost_gradient_from_reference();
    make_const_cost_from_reference();

    // LTV update
    if (LTV)
    {
        update_equality_constraints_LTV();
        solver_qp.set_C_vec(C_i_vec);
    }
        
    // update solver
    solver_qp.set_crhs_vec(crhs_i_vec);
    solver_qp.set_q_vec(q_i_vec);
    solver_qp.set_b(b);

    // solve
    bool feasible = solve_optimization_problem();

    // return control input for time zero (or 1 if using u1 control)
    if (u1_control)
        return std::make_pair(u_vec[1], feasible);
    else
        return std::make_pair(u_vec[0], feasible);
}

// get methods
std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>> MPC_QP::get_trajectory()
{
    return std::make_pair(x_vec, u_vec);
}

double MPC_QP::get_objective()
{
    return objective;
}

// get problem dimensions
void MPC_QP::get_problem_dimensions()
{
    // full problem dimensions
    n_y_tot = 0; // init
    n_ineq_tot = 0; // init
    for (int k=0; k<=n_horizon; ++k)
    {
        if (q_i_vec.count(k))
            n_y_vec[k] = q_i_vec[k].size();
        else
            n_y_vec[k] = q_i_nom.size();
        n_y_tot += n_y_vec[k];

        if (w_i_vec.count(k))
            n_ineq_vec[k] = w_i_vec[k].size();
        else
            n_ineq_vec[k] = w_i_nom.size();
        n_ineq_tot += n_ineq_vec[k];
    }

    n_eq_tot = 0; // init
    for (int k=0; k<n_horizon; ++k)
    {
        if (crhs_i_vec.count(k))
            n_eq_vec[k] = crhs_i_vec[k].size();
        else
            n_eq_vec[k] = crhs_i_nom.size();
        n_eq_tot += n_eq_vec[k];
    }

    // indexing
    int idx_offset = 0; // init
    int idx_eq_offset = 0; // init
    int idx_ineq_offset = 0; // init
    for (int k=0; k<=n_horizon; ++k)
    {
        idx_state[k] = get_indices(idx_offset, n_state);
        idx_input[k] = get_indices(idx_offset + n_state, n_input);
        idx_y[k] = get_indices(idx_offset, n_y_vec[k]);
        idx_ineq[k] = get_indices(idx_ineq_offset, n_ineq_vec[k]);
        idx_eq[k] = get_indices(idx_eq_offset, n_eq_vec[k]);

        idx_offset += n_y_vec[k];
        idx_ineq_offset += n_ineq_vec[k];
        idx_eq_offset += n_eq_vec[k];
    }
}

// configure solver
void MPC_QP::configure_solver()
{
    // setup
    solver_qp.setup(P_i_nom, q_i_nom, C_i_nom, D_i_nom, crhs_i_nom, G_i_nom, w_i_nom, n_horizon, b);
    solver_qp.set_P_vec(P_i_vec);
    solver_qp.set_q_vec(q_i_vec);
    solver_qp.set_C_vec(C_i_vec);
    solver_qp.set_D_vec(D_i_vec);
    solver_qp.set_crhs_vec(crhs_i_vec);
    solver_qp.set_G_vec(G_i_vec);
    solver_qp.set_w_vec(w_i_vec);

    // set solver settings
    solver_qp.set_settings(qp_settings);
}

// solve optimization problem
bool MPC_QP::solve_optimization_problem()
{
    // solve
    QP_IP_MPC::QP_results results = solver_qp.solve();
    solution = results.y;
    objective = results.objective;

    // check if solution is valid
    bool feasible = (results.converged && results.feas);
    
    // get state and input trajectories
    get_state_input_trajectories();

    return feasible;
}

// get state and input trajectories
void MPC_QP::get_state_input_trajectories()
{
    // clear vectors
    x_vec.clear();
    u_vec.clear();

    // get state and input trajectories
    for (int k=0; k<=n_horizon; ++k)
    {
        x_vec.push_back(solution.segment(idx_state[k][0], n_state));
        if (k < n_horizon)
            u_vec.push_back(solution.segment(idx_input[k][0], n_input));
    }
}

// build problem matrices

// make inequality constraints
void MPC_QP::make_inequality_constraints()
{
    // declarations
    int m_Gx, n_Gx, m_Gu, n_Gu;
    std::vector<Eigen::Triplet<double>> tripvec_G_i_nom, tripvec_G_term;

    if (Hrep_constraints)
    {
        // matrix dimensions
        m_Gx = Ax_ineq.rows();
        n_Gx = Ax_ineq.cols();
        m_Gu = Au_ineq.rows();
        n_Gu = Au_ineq.cols();

        // build nominal inequality constraints
        // state vector: y = [x, u, sx_ineq, su_ineq]

        tripvec_G_i_nom.reserve(Ax_ineq.nonZeros() + Au_ineq.nonZeros() + 2*m_Gx + 2*m_Gu);

        if (soft_state_ineq && soft_input_ineq)
        {
            get_triplets_offset(Ax_ineq, tripvec_G_i_nom, 0, 0);
            get_triplets_offset_neg_identity(m_Gx, tripvec_G_i_nom, 0, n_Gx + n_Gu);
            get_triplets_offset(Au_ineq, tripvec_G_i_nom, m_Gx, n_Gx);
            get_triplets_offset_neg_identity(m_Gu, tripvec_G_i_nom, m_Gx, n_Gx + n_Gu + n_Gx);
            get_triplets_offset_identity(m_Gx, tripvec_G_i_nom, m_Gx + m_Gu, n_Gx + n_Gu);
            get_triplets_offset_neg_identity(m_Gu, tripvec_G_i_nom, m_Gx + m_Gu + m_Gx, n_Gx + n_Gu);
            get_triplets_offset_identity(m_Gu, tripvec_G_i_nom, m_Gx + m_Gu + 2*m_Gx, n_Gx + n_Gu + m_Gx);
            get_triplets_offset_neg_identity(m_Gu, tripvec_G_i_nom, m_Gx + m_Gu + 2*m_Gx + m_Gu, n_Gx + n_Gu + m_Gx);

            G_i_nom.resize(m_Gx + m_Gu + 2*m_Gx + 2*m_Gu, n_Gx + n_Gu + m_Gx + m_Gu);
            G_i_nom.setFromTriplets(tripvec_G_i_nom.begin(), tripvec_G_i_nom.end());

            w_i_nom = Eigen::VectorXd::Zero(m_Gx + m_Gu + 2*m_Gx + 2*m_Gu);
            w_i_nom.segment(0, m_Gx) = bx_ineq;
            w_i_nom.segment(m_Gx, m_Gu) = bu_ineq;
            w_i_nom.segment(m_Gx + m_Gu, m_Gx) = sigma_max_x;
            w_i_nom.segment(m_Gx + m_Gu + m_Gx, m_Gx) = Eigen::VectorXd::Zero(sigma_max_x.size());
            w_i_nom.segment(m_Gx + m_Gu + 2*m_Gx, m_Gu) = sigma_max_u;
            w_i_nom.segment(m_Gx + m_Gu + 2*m_Gx + m_Gu, m_Gu) = Eigen::VectorXd::Zero(sigma_max_u.size());
        }
        else if (soft_state_ineq)
        {
            get_triplets_offset(Ax_ineq, tripvec_G_i_nom, 0, 0);
            get_triplets_offset_neg_identity(m_Gx, tripvec_G_i_nom, 0, n_Gx + n_Gu);
            get_triplets_offset(Au_ineq, tripvec_G_i_nom, m_Gx, n_Gx);
            get_triplets_offset_identity(m_Gx, tripvec_G_i_nom, m_Gx + m_Gu, n_Gx + n_Gu);
            get_triplets_offset_neg_identity(m_Gx, tripvec_G_i_nom, m_Gx + m_Gu + m_Gx, n_Gx + n_Gu);

            G_i_nom.resize(m_Gx + m_Gu + 2*m_Gx, n_Gx + n_Gu + m_Gx);
            G_i_nom.setFromTriplets(tripvec_G_i_nom.begin(), tripvec_G_i_nom.end());

            w_i_nom = Eigen::VectorXd::Zero(m_Gx + m_Gu + 2*m_Gx);
            w_i_nom.segment(0, m_Gx) = bx_ineq;
            w_i_nom.segment(m_Gx, m_Gu) = bu_ineq;
            w_i_nom.segment(m_Gx + m_Gu, m_Gx) = sigma_max_x;
            w_i_nom.segment(m_Gx + m_Gu + m_Gx, m_Gx) = Eigen::VectorXd::Zero(sigma_max_x.size());        
        }
        else if (soft_input_ineq)
        {
            get_triplets_offset(Ax_ineq, tripvec_G_i_nom, 0, 0);
            get_triplets_offset(Au_ineq, tripvec_G_i_nom, m_Gx, n_Gx);
            get_triplets_offset_neg_identity(m_Gu, tripvec_G_i_nom, m_Gx, n_Gx + n_Gu);
            get_triplets_offset_identity(m_Gu, tripvec_G_i_nom, m_Gx + m_Gu, n_Gx + n_Gu);
            get_triplets_offset_neg_identity(m_Gu, tripvec_G_i_nom, m_Gx + m_Gu + m_Gu, n_Gx + n_Gu);

            G_i_nom.resize(m_Gx + m_Gu + 2*m_Gu, n_Gx + n_Gu + m_Gu);
            G_i_nom.setFromTriplets(tripvec_G_i_nom.begin(), tripvec_G_i_nom.end());

            w_i_nom = Eigen::VectorXd::Zero(m_Gx + m_Gu + 2*m_Gu);
            w_i_nom.segment(0, m_Gx) = bx_ineq;
            w_i_nom.segment(m_Gx, m_Gu) = bu_ineq;
            w_i_nom.segment(m_Gx + m_Gu, m_Gu) = sigma_max_u;
            w_i_nom.segment(m_Gx + m_Gu + m_Gu, m_Gu) = Eigen::VectorXd::Zero(sigma_max_u.size());
        }
        else
        {
            get_triplets_offset(Ax_ineq, tripvec_G_i_nom, 0, 0);
            get_triplets_offset(Au_ineq, tripvec_G_i_nom, m_Gx, n_Gx);

            G_i_nom.resize(m_Gx + m_Gu, n_Gx + n_Gu);
            G_i_nom.setFromTriplets(tripvec_G_i_nom.begin(), tripvec_G_i_nom.end());

            w_i_nom = Eigen::VectorXd::Zero(m_Gx + m_Gu);
            w_i_nom.segment(0, m_Gx) = bx_ineq;
            w_i_nom.segment(m_Gx, m_Gu) = bu_ineq;
        }

        // reset off-nominal inequality constraints
        G_i_vec.clear();
        w_i_vec.clear();

        // build terminal constraints if applicable
        if (terminal_state_ineq_specified)
        {
            m_Gx = Ax_term_ineq.rows();
            n_Gx = Ax_term_ineq.cols();
            tripvec_G_term.reserve(Ax_term_ineq.nonZeros() + 2*m_Gx);

            if (soft_term_state_ineq)
            {
                get_triplets_offset(Ax_term_ineq, tripvec_G_term, 0, 0);
                get_triplets_offset_neg_identity(m_Gx, tripvec_G_term, 0, n_Gx);
                get_triplets_offset_identity(m_Gx, tripvec_G_term, m_Gx, n_Gx);
                get_triplets_offset_neg_identity(m_Gx, tripvec_G_term, m_Gx + m_Gx, n_Gx);

                G_i_vec[n_horizon].resize(m_Gx + 2*m_Gx, n_Gx + m_Gx);
                G_i_vec[n_horizon].setFromTriplets(tripvec_G_term.begin(), tripvec_G_term.end());

                w_i_vec[n_horizon] = Eigen::VectorXd::Zero(m_Gx + 2*m_Gx);
                w_i_vec[n_horizon].segment(0, m_Gx) = bx_term_ineq;
                w_i_vec[n_horizon].segment(m_Gx, m_Gx) = sigma_max_x_term;
                w_i_vec[n_horizon].segment(m_Gx + m_Gx, m_Gx) = Eigen::VectorXd::Zero(sigma_max_x_term.size());
            }
            else
            {
                get_triplets_offset(Ax_term_ineq, tripvec_G_term, 0, 0);

                G_i_vec[n_horizon].resize(m_Gx, n_Gx);
                G_i_vec[n_horizon].setFromTriplets(tripvec_G_term.begin(), tripvec_G_term.end());

                w_i_vec[n_horizon] = bx_term_ineq;
            }
        }
        else
        {
            m_Gx = Ax_ineq.rows();
            n_Gx = Ax_ineq.cols();
            tripvec_G_term.reserve(Ax_ineq.nonZeros() + 2*m_Gx);

            if (soft_state_ineq)
            {
                get_triplets_offset(Ax_ineq, tripvec_G_term, 0, 0);
                get_triplets_offset_neg_identity(m_Gx, tripvec_G_term, 0, n_Gx);
                get_triplets_offset_identity(m_Gx, tripvec_G_term, m_Gx, n_Gx);
                get_triplets_offset_neg_identity(m_Gx, tripvec_G_term, m_Gx + m_Gx, n_Gx);

                G_i_vec[n_horizon].resize(m_Gx + 2*m_Gx, n_Gx + m_Gx);
                G_i_vec[n_horizon].setFromTriplets(tripvec_G_term.begin(), tripvec_G_term.end());

                w_i_vec[n_horizon] = Eigen::VectorXd::Zero(m_Gx + 2*m_Gx);
                w_i_vec[n_horizon].segment(0, m_Gx) = bx_ineq;
                w_i_vec[n_horizon].segment(m_Gx, m_Gx) = sigma_max_x;
                w_i_vec[n_horizon].segment(m_Gx + m_Gx, m_Gx) = Eigen::VectorXd::Zero(sigma_max_x.size());
            }
            else
            {
                get_triplets_offset(Ax_ineq, tripvec_G_term, 0, 0);

                G_i_vec[n_horizon].resize(m_Gx, n_Gx);
                G_i_vec[n_horizon].setFromTriplets(tripvec_G_term.begin(), tripvec_G_term.end());

                w_i_vec[n_horizon] = bx_ineq;
            }
        }
    } // end if (Hrep_constraints)
    else // conzono constraints
    {
        // nominal inequality constraints

        // dimensions
        int nx_soft = soft_state_ineq ? Zc_x.n : 0;
        int nu_soft = soft_input_ineq ? Zc_u.n : 0;

        int n_box_cons = 2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG + 2*nx_soft + 2*nu_soft;
        tripvec_G_i_nom.reserve(n_box_cons);
        w_i_nom = Eigen::VectorXd::Zero(n_box_cons);

        // hard state constraints
        get_triplets_offset_identity(n_state, tripvec_G_i_nom, 0, 0);
        get_triplets_offset_neg_identity(n_state, tripvec_G_i_nom, n_state, 0);
        w_i_nom.segment(0, n_state) = x_max_hard;
        w_i_nom.segment(n_state, n_state) = -x_min_hard;

        // hard input constraints
        get_triplets_offset_identity(n_input, tripvec_G_i_nom, 2*n_state, n_state);
        get_triplets_offset_neg_identity(n_input, tripvec_G_i_nom, 2*n_state+n_input, n_state);
        w_i_nom.segment(2*n_state, n_input) = u_max_hard;
        w_i_nom.segment(2*n_state+n_input, n_input) = -u_min_hard;

        // conzono box constraints
        get_triplets_offset_identity(Zc_x.nG, tripvec_G_i_nom, 2*n_state + 2*n_input, n_state + n_input);
        get_triplets_offset_neg_identity(Zc_x.nG, tripvec_G_i_nom, 2*n_state + 2*n_input + Zc_x.nG, n_state + n_input);
        w_i_nom.segment(2*n_state + 2*n_input, Zc_x.nG) = Eigen::VectorXd::Ones(Zc_x.nG);
        w_i_nom.segment(2*n_state + 2*n_input + Zc_x.nG, Zc_x.nG) = Eigen::VectorXd::Zero(Zc_x.nG);

        get_triplets_offset_identity(Zc_u.nG, tripvec_G_i_nom, 2*n_state + 2*n_input + 2*Zc_x.nG, n_state + n_input + Zc_x.nG);
        get_triplets_offset_neg_identity(Zc_u.nG, tripvec_G_i_nom, 2*n_state + 2*n_input + 2*Zc_x.nG + Zc_u.nG, n_state + n_input + Zc_x.nG);
        w_i_nom.segment(2*n_state + 2*n_input + 2*Zc_x.nG, Zc_u.nG) = Eigen::VectorXd::Ones(Zc_u.nG);
        w_i_nom.segment(2*n_state + 2*n_input + 2*Zc_x.nG + Zc_u.nG, Zc_u.nG) = Eigen::VectorXd::Zero(Zc_u.nG);

        // slack vars
        get_triplets_offset_identity(nx_soft, tripvec_G_i_nom, 2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG, 
            n_state + n_input + Zc_x.nG + Zc_u.nG);
        get_triplets_offset_neg_identity(nx_soft, tripvec_G_i_nom, 2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG + nx_soft,
            n_state + n_input + Zc_x.nG + Zc_u.nG);
        w_i_nom.segment(2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG, nx_soft) = sigma_max_x;
        w_i_nom.segment(2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG + nx_soft, nx_soft) = sigma_max_x;

        get_triplets_offset_identity(nu_soft, tripvec_G_i_nom, 2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG + 2*nx_soft, 
            n_state + n_input + Zc_x.nG + Zc_u.nG + nx_soft);
        get_triplets_offset_neg_identity(nu_soft, tripvec_G_i_nom, 2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG + 2*nx_soft + nu_soft,
            n_state + n_input + Zc_x.nG + Zc_u.nG + nx_soft);
        w_i_nom.segment(2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG + 2*nx_soft, nu_soft) = sigma_max_u;
        w_i_nom.segment(2*n_state + 2*n_input + 2*Zc_x.nG + 2*Zc_u.nG + 2*nx_soft + nu_soft, nu_soft) = sigma_max_u;

        // build constraint matrix
        G_i_nom.resize(n_box_cons, n_box_cons/2);
        G_i_nom.setFromTriplets(tripvec_G_i_nom.begin(), tripvec_G_i_nom.end());

        // reset off-nominal inequality constraints
        G_i_vec.clear();
        w_i_vec.clear();

        // build terminal constraints if applicable
        ZonoCpp::ConZono * Zc_x_term_ptr = nullptr;
        Eigen::VectorXd * x_max_ptr = nullptr;
        Eigen::VectorXd * x_min_ptr = nullptr;
        Eigen::VectorXd * sigma_max_x_ptr = nullptr;
        if (terminal_state_ineq_specified)
        {
            Zc_x_term_ptr = &Zc_x_term;
            x_max_ptr = &x_term_max_hard;
            x_min_ptr = &x_term_min_hard;
            sigma_max_x_ptr = &sigma_max_x_term;
        }
        else
        {
            Zc_x_term_ptr = &Zc_x;
            x_max_ptr = &x_max_hard;
            x_min_ptr = &x_min_hard;
            sigma_max_x_ptr = &sigma_max_x;
        }

        int nx_term_soft = soft_term_state_ineq ? Zc_x_term_ptr->n : 0;
        n_box_cons = 2*n_state + 2*Zc_x_term_ptr->nG + 2*nx_term_soft;
        
        tripvec_G_term.reserve(n_box_cons);
        w_i_vec[n_horizon] = Eigen::VectorXd::Zero(n_box_cons);

        // hard state constraints
        get_triplets_offset_identity(n_state, tripvec_G_term, 0, 0);
        get_triplets_offset_neg_identity(n_state, tripvec_G_term, n_state, 0);
        w_i_vec[n_horizon].segment(0, n_state) = *x_max_ptr;
        w_i_vec[n_horizon].segment(n_state, n_state) = -1.0*(*x_min_ptr);

        // conzono box constraints
        get_triplets_offset_identity(Zc_x_term_ptr->nG, tripvec_G_term, 2*n_state, n_state);
        get_triplets_offset_neg_identity(Zc_x_term_ptr->nG, tripvec_G_term, 2*n_state + Zc_x_term_ptr->nG, n_state);
        w_i_vec[n_horizon].segment(2*n_state, Zc_x_term_ptr->nG) = Eigen::VectorXd::Ones(Zc_x_term_ptr->nG);
        w_i_vec[n_horizon].segment(2*n_state + Zc_x_term_ptr->nG, Zc_x_term_ptr->nG) = Eigen::VectorXd::Zero(Zc_x_term_ptr->nG);

        // slack vars
        get_triplets_offset_identity(nx_term_soft, tripvec_G_term, 2*n_state + 2*Zc_x_term_ptr->nG, n_state + Zc_x_term_ptr->nG);
        get_triplets_offset_neg_identity(nx_term_soft, tripvec_G_term, 2*n_state + 2*Zc_x_term_ptr->nG + nx_term_soft, n_state + Zc_x_term_ptr->nG);
        w_i_vec[n_horizon].segment(2*n_state + 2*Zc_x_term_ptr->nG, nx_term_soft) = *sigma_max_x_ptr;
        w_i_vec[n_horizon].segment(2*n_state + 2*Zc_x_term_ptr->nG + nx_term_soft, nx_term_soft) = *sigma_max_x_ptr;

        // build matrix
        G_i_vec[n_horizon].resize(n_box_cons, n_box_cons/2);
        G_i_vec[n_horizon].setFromTriplets(tripvec_G_term.begin(), tripvec_G_term.end());
    }
}

// make equality constraints
void MPC_QP::make_equality_constraints()
{
    // state vector y = [x, u, sx_ineq, su_ineq]

    // declarations
    int m_Gx, n_Gx, m_Gu, n_Gu, m_A, n_A, m_B, n_B, n_cons_init;
    std::vector<Eigen::Triplet<double>> tripvec_C_i_nom, tripvec_D_i_nom, tripvec_C_0, tripvec_D_1, tripvec_D_N;

    // H-rep constraint case
    if (Hrep_constraints)
    {
        // make nominal equality matrices / vector
        // matrix dimensions
        m_Gx = Ax_ineq.rows();
        n_Gx = Ax_ineq.cols();
        m_Gu = Au_ineq.rows();
        n_Gu = Au_ineq.cols();
        m_A = A_dyn.rows();
        n_A = A_dyn.cols();
        m_B = B_dyn.rows();
        n_B = B_dyn.cols();

        // nominal case
        tripvec_C_i_nom.reserve(A_dyn.nonZeros() + B_dyn.nonZeros());
        tripvec_D_i_nom.reserve(m_A);

        get_triplets_offset(A_dyn, tripvec_C_i_nom, 0, 0);
        get_triplets_offset(B_dyn, tripvec_C_i_nom, 0, n_A);
        get_triplets_offset_neg_identity(m_A, tripvec_D_i_nom, 0, 0);

        if (soft_state_ineq && soft_input_ineq)
        {
            C_i_nom.resize(m_A, n_A + n_B + m_Gx + m_Gu);
            D_i_nom.resize(m_A, n_A + n_B + m_Gx + m_Gu);
        }
        else if (soft_state_ineq)
        {
            C_i_nom.resize(m_A, n_A + n_B + m_Gx);
            D_i_nom.resize(m_A, n_A + n_B + m_Gx);
        }
        else if (soft_input_ineq)
        {
            C_i_nom.resize(m_A, n_A + n_B + m_Gu);
            D_i_nom.resize(m_A, n_A + n_B + m_Gu);
        }
        else
        {
            C_i_nom.resize(m_A, n_A + n_B);
            D_i_nom.resize(m_A, n_A + n_B);
        }

        C_i_nom.setFromTriplets(tripvec_C_i_nom.begin(), tripvec_C_i_nom.end());
        D_i_nom.setFromTriplets(tripvec_D_i_nom.begin(), tripvec_D_i_nom.end());

        crhs_i_nom = Eigen::VectorXd::Zero(m_A);

        // indexing for dynamics equation
        idx_dyn_nom = get_indices(0, m_A);

        // reset off-nominal constraints
        C_i_vec.clear();
        D_i_vec.clear();
        crhs_i_vec.clear();

        // initial constraints
        tripvec_C_0.reserve(A_dyn.nonZeros() + B_dyn.nonZeros() + m_A + n_B);
        tripvec_D_1.reserve(m_A);

        if (u1_control)
        {
            get_triplets_offset_neg_identity(m_A, tripvec_C_0, 0, 0);
            get_triplets_offset_neg_identity(n_B, tripvec_C_0, m_A, m_A);
            get_triplets_offset(A_dyn, tripvec_C_0, m_A + n_B, 0);
            get_triplets_offset(B_dyn, tripvec_C_0, m_A + n_B, n_A);
            get_triplets_offset_neg_identity(m_A, tripvec_D_1, m_A + n_B, 0);
        }
        else
        {
            get_triplets_offset_neg_identity(m_A, tripvec_C_0, 0, 0);
            get_triplets_offset(A_dyn, tripvec_C_0, m_A, 0);
            get_triplets_offset(B_dyn, tripvec_C_0, m_A, n_A);
            get_triplets_offset_neg_identity(m_A, tripvec_D_1, m_A, 0);
        }

        n_cons_init = u1_control ? 2 * m_A + n_B : 2 * m_A;
        if (soft_state_ineq && soft_input_ineq)
        {
            C_i_vec[0].resize(n_cons_init, n_A + n_B + m_Gx + m_Gu);
            D_i_vec[1].resize(n_cons_init, n_A + n_B + m_Gx + m_Gu);
        }
        else if (soft_state_ineq)
        {
            C_i_vec[0].resize(n_cons_init, n_A + n_B + m_Gx);
            D_i_vec[1].resize(n_cons_init, n_A + n_B + m_Gx);
        }
        else if (soft_input_ineq)
        {
            C_i_vec[0].resize(n_cons_init, n_A + n_B + m_Gu);
            D_i_vec[1].resize(n_cons_init, n_A + n_B + m_Gu);
        }
        else
        {
            C_i_vec[0].resize(n_cons_init, n_A + n_B);
            D_i_vec[1].resize(n_cons_init, n_A + n_B);
        }

        C_i_vec[0].setFromTriplets(tripvec_C_0.begin(), tripvec_C_0.end());
        D_i_vec[1].setFromTriplets(tripvec_D_1.begin(), tripvec_D_1.end());

        crhs_i_vec[0] = Eigen::VectorXd::Zero(C_i_vec[0].rows());

        // indexing for dynamics
        if (u1_control)
            idx_dyn_vec[0] = get_indices(m_A+n_B, m_A);
        else
            idx_dyn_vec[0] = get_indices(m_A, m_A);

        // terminal constraints
        tripvec_D_N.reserve(m_A);

        get_triplets_offset_neg_identity(m_A, tripvec_D_N, 0, 0);

        if (terminal_state_ineq_specified)
            m_Gx = Ax_term_ineq.rows();
        else
            m_Gx = Ax_ineq.rows();

        if (soft_term_state_ineq)
            D_i_vec[n_horizon].resize(m_A, m_A + m_Gx);
        else
            D_i_vec[n_horizon].resize(m_A, m_A);

        D_i_vec[n_horizon].setFromTriplets(tripvec_D_N.begin(), tripvec_D_N.end());

    } // end if (Hrep_constraints)
    else // conzono constraints
    {
        // make nominal equality matrices / vector
        // matrix dimensions
        m_A = A_dyn.rows();
        n_A = A_dyn.cols();
        m_B = B_dyn.rows();
        n_B = B_dyn.cols();

        // indexing for dynamics equation
        idx_dyn_nom = get_indices(0, m_A);

        // nominal case
        tripvec_C_i_nom.reserve(A_dyn.nonZeros() + B_dyn.nonZeros() + Zc_x.G.nonZeros() + Zc_u.G.nonZeros() + Zc_x.n + Zc_u.n);
        tripvec_D_i_nom.reserve(m_A);

        // dynamics
        get_triplets_offset(A_dyn, tripvec_C_i_nom, 0, 0);
        get_triplets_offset(B_dyn, tripvec_C_i_nom, 0, n_A);
        get_triplets_offset_neg_identity(m_A, tripvec_D_i_nom, 0, 0);

        // state/input constraints
        get_triplets_offset(Pp_x, tripvec_C_i_nom, m_A, 0);
        get_triplets_offset(-1*Zc_x.G, tripvec_C_i_nom, m_A, n_A + n_B);
        get_triplets_offset(Pp_u, tripvec_C_i_nom, m_A + Zc_x.n, n_A);
        get_triplets_offset(-1*Zc_u.G, tripvec_C_i_nom, m_A + Zc_x.n, n_A + n_B + Zc_x.nG);

        // softening
        int nx_soft = soft_state_ineq ? Zc_x.n : 0;
        int nu_soft = soft_input_ineq ? Zc_u.n : 0;

        get_triplets_offset_identity(nx_soft, tripvec_C_i_nom, m_A, n_A + n_B + Zc_x.nG + Zc_u.nG); // state constraint softening
        get_triplets_offset_identity(nu_soft, tripvec_C_i_nom, m_A + Zc_x.n, n_A + n_B + Zc_x.nG + Zc_u.nG + nx_soft); // input constraint softening

        // conzono constraints
        get_triplets_offset(Zc_x.A, tripvec_C_i_nom, m_A + Zc_x.n + Zc_u.n, n_A + n_B);
        get_triplets_offset(Zc_u.A, tripvec_C_i_nom, m_A + Zc_x.n + Zc_u.n + Zc_x.nC, n_A + n_B + Zc_x.nG);

        // build nominal constraints
        C_i_nom.resize(m_A + Zc_x.n + Zc_u.n + Zc_x.nC + Zc_u.nC, n_A + n_B + Zc_x.nG + Zc_u.nG + nx_soft + nu_soft);
        D_i_nom.resize(m_A + Zc_x.n + Zc_u.n + Zc_x.nC + Zc_u.nC, n_A + n_B + Zc_x.nG + Zc_u.nG + nx_soft + nu_soft);

        C_i_nom.setFromTriplets(tripvec_C_i_nom.begin(), tripvec_C_i_nom.end());
        D_i_nom.setFromTriplets(tripvec_D_i_nom.begin(), tripvec_D_i_nom.end());

        crhs_i_nom = Eigen::VectorXd::Zero(C_i_nom.rows());
        crhs_i_nom.segment(m_A, Zc_x.n) = -1*Zc_x.c;
        crhs_i_nom.segment(m_A+Zc_x.n, Zc_u.n) = -1*Zc_u.c;
        crhs_i_nom.segment(m_A+Zc_x.n+Zc_u.n, Zc_x.nC) = -1*Zc_x.b;
        crhs_i_nom.segment(m_A+Zc_x.n+Zc_u.n+Zc_x.nC, Zc_u.nC) = -1*Zc_u.b;

        // reset off-nominal constraints
        C_i_vec.clear();
        D_i_vec.clear();
        crhs_i_vec.clear();

        // initial constraints
        tripvec_C_0.reserve(C_i_nom.nonZeros() + m_A + n_B);
        tripvec_D_1.reserve(D_i_nom.nonZeros());

        if (u1_control)
        {
            // initial state and input
            get_triplets_offset_neg_identity(m_A, tripvec_C_0, 0, 0);
            get_triplets_offset_neg_identity(n_B, tripvec_C_0, m_A, m_A);
            
            // nominal equality constraints
            get_triplets_offset(C_i_nom, tripvec_C_0, m_A + n_B, 0);
            get_triplets_offset(D_i_nom, tripvec_D_1, m_A + n_B, 0);
        }
        else
        {
            // initial state
            get_triplets_offset_neg_identity(m_A, tripvec_C_0, 0, 0);
            
            // nominal equality constraints
            get_triplets_offset(C_i_nom, tripvec_C_0, m_A, 0);
            get_triplets_offset(D_i_nom, tripvec_D_1, m_A, 0);
        }

        n_cons_init = (u1_control ? m_A + n_B : m_A) + C_i_nom.rows();
        C_i_vec[0].resize(n_cons_init, C_i_nom.cols());
        D_i_vec[1].resize(n_cons_init, D_i_nom.cols());

        C_i_vec[0].setFromTriplets(tripvec_C_0.begin(), tripvec_C_0.end());
        D_i_vec[1].setFromTriplets(tripvec_D_1.begin(), tripvec_D_1.end());
        
        int n_new_cons = u1_control ? m_A + n_B : m_A;
        crhs_i_vec[0] = Eigen::VectorXd::Zero(n_new_cons + C_i_nom.rows());
        crhs_i_vec[0].segment(n_new_cons, C_i_nom.rows()) = crhs_i_nom;

        // indexing for dynamics
        if (u1_control)
            idx_dyn_vec[0] = get_indices(m_A+n_B, m_A);
        else
            idx_dyn_vec[0] = get_indices(m_A, m_A);

        // terminal constraints
        ZonoCpp::ConZono * Zc_x_term_ptr = nullptr;
        Eigen::SparseMatrix<double> * Pp_x_term_ptr = nullptr;
        if (terminal_state_ineq_specified)
        {
            Zc_x_term_ptr = &Zc_x_term;
            Pp_x_term_ptr = &Pp_x_term;
        }
        else
        {    
            Zc_x_term_ptr = &Zc_x;
            Pp_x_term_ptr = &Pp_x;
        }

        tripvec_D_N.reserve(D_i_nom.nonZeros() + Zc_x_term_ptr->G.nonZeros());

        get_triplets_offset_neg_identity(m_A, tripvec_D_N, 0, 0);
        get_triplets_offset(*Pp_x_term_ptr, tripvec_D_N, D_i_nom.rows(), 0);
        get_triplets_offset(-1*Zc_x_term_ptr->G, tripvec_D_N, D_i_nom.rows(), n_A); // no input constraints

        int nx_term_soft = soft_term_state_ineq ? Zc_x_term_ptr->n : 0;
        get_triplets_offset_identity(nx_term_soft, tripvec_D_N, D_i_nom.rows(), n_A + Zc_x_term_ptr->nG);

        get_triplets_offset(Zc_x_term_ptr->A, tripvec_D_N, D_i_nom.rows() + Zc_x_term_ptr->n, n_A);

        D_i_vec[n_horizon].resize(D_i_nom.rows() + Zc_x_term_ptr->n + Zc_x_term_ptr->nC, n_A + Zc_x_term_ptr->nG + nx_term_soft);
        D_i_vec[n_horizon].setFromTriplets(tripvec_D_N.begin(), tripvec_D_N.end());

        if (!C_i_vec.count(n_horizon-1))
            C_i_vec[n_horizon-1] = C_i_nom;
        C_i_vec[n_horizon-1].conservativeResize(D_i_vec[n_horizon].rows(), C_i_vec[n_horizon-1].cols());

        if (!crhs_i_vec.count(n_horizon-1))
            crhs_i_vec[n_horizon-1] = crhs_i_nom;
        crhs_i_vec[n_horizon-1].conservativeResize(D_i_vec[n_horizon].rows());
        crhs_i_vec[n_horizon-1].segment(D_i_nom.rows(), Zc_x_term_ptr->n) = -1*Zc_x_term_ptr->c;
        crhs_i_vec[n_horizon-1].segment(D_i_nom.rows() + Zc_x_term_ptr->n, Zc_x_term_ptr->nC) = -1*Zc_x_term_ptr->b;

        // track indexing
        start_idx_term_cons = D_i_nom.rows();
    }

    // if LTV, need to do constraints by time step
    if (LTV)
    {
        // loop through time steps
        for (int k = 1; k < n_horizon; ++k)
        {
            if (!C_i_vec.count(k))
                C_i_vec[k] = C_i_nom;
            update_sparse_matrix(C_i_vec[k], A_dyn_vec[k], 0, 0);
            update_sparse_matrix(C_i_vec[k], B_dyn_vec[k], 0, n_state);
        }
    }

    // make IC constraints
    make_IC_constraint();
}

// LTV update
void MPC_QP::update_equality_constraints_LTV()
{
    // initial time step
    int i_offset_0 = u1_control ? n_state + n_input : n_state;
    if (!C_i_vec.count(0))
        C_i_vec[0] = C_i_nom;

    update_sparse_matrix(C_i_vec[0], A_dyn_vec[0], i_offset_0, 0);
    update_sparse_matrix(C_i_vec[0], B_dyn_vec[0], i_offset_0, n_state);

    // loop through other time steps
    for (int k=1; k < n_horizon; ++k)
    {
        if (!C_i_vec.count(k))
            C_i_vec[k] = C_i_nom;
        update_sparse_matrix(C_i_vec[k], A_dyn_vec[k], 0, 0);
        update_sparse_matrix(C_i_vec[k], B_dyn_vec[k], 0, n_state);
    }
}

// terminal constraint update
void MPC_QP::update_ref_dep_term_constraint()
{
    // only support conzono constraints at the moment
    if (!conzono_constraints)
        return;

    // get conzono terminal constraint
    ZonoCpp::ConZono * Zc_x_term_ptr = terminal_state_ineq_specified ? &Zc_x_term : &Zc_x;

    // Minkowski sum with state reference
    crhs_i_vec[n_horizon-1].segment(start_idx_term_cons, n_state) = -(Zc_x_term_ptr->c + x_ref.back());
}

// make IC constraint
void MPC_QP::make_IC_constraint()
{
    // get original crhs at step 0
    if (!crhs_i_vec.count(0))
        crhs_i_vec[0] = crhs_i_nom;

    // update initial condition
    crhs_i_vec[0].segment(0, n_state) = x0;

    // if using u1 control, update u0 constraint
    if (u1_control)
        crhs_i_vec[0].segment(n_state, n_input) = u0;
}

// make cost matrices
void MPC_QP::make_cost_function()
{
    // declare
    std::vector<Eigen::Triplet<double>> tripvec_P_i_nom, tripvec_P_N;
    int offset;
    Eigen::SparseMatrix<double> * P = nullptr;
    if (terminal_cost_specified)
        P = &P_cost;
    else
        P = &Q_cost;

    if (Hrep_constraints)
    {
        // state vector y = [x, u, sx_ineq, su_ineq]

        // make nominal cost matrices / vector

        // nominal case
        tripvec_P_i_nom.reserve(Q_cost.nonZeros() + R_cost.nonZeros() + 2*N_cost.nonZeros() 
            + Qx_ineq_cost.nonZeros() + Qu_ineq_cost.nonZeros());

        offset = 0;
        if (soft_state_ineq && soft_input_ineq)
        {
            get_triplets_offset(Q_cost, tripvec_P_i_nom, 0, 0);
            get_triplets_offset(N_cost, tripvec_P_i_nom, 0, Q_cost.cols());
            offset += Q_cost.rows();

            get_triplets_offset(R_cost, tripvec_P_i_nom, offset, offset);
            get_triplets_offset(N_cost.transpose(), tripvec_P_i_nom, offset, 0);
            offset += R_cost.rows();

            get_triplets_offset(Qx_ineq_cost, tripvec_P_i_nom, offset, offset);
            offset += Qx_ineq_cost.rows();

            get_triplets_offset(Qu_ineq_cost, tripvec_P_i_nom, offset, offset);
            offset += Qu_ineq_cost.rows();

            P_i_nom.resize(offset, offset);
        }
        else if (soft_state_ineq)
        {
            get_triplets_offset(Q_cost, tripvec_P_i_nom, 0, 0);
            get_triplets_offset(N_cost, tripvec_P_i_nom, 0, Q_cost.cols());
            offset += Q_cost.rows();

            get_triplets_offset(R_cost, tripvec_P_i_nom, offset, offset);
            get_triplets_offset(N_cost.transpose(), tripvec_P_i_nom, offset, 0);
            offset += R_cost.rows();

            get_triplets_offset(Qx_ineq_cost, tripvec_P_i_nom, offset, offset);
            offset += Qx_ineq_cost.rows();

            P_i_nom.resize(offset, offset);
        }
        else if (soft_input_ineq)
        {
            get_triplets_offset(Q_cost, tripvec_P_i_nom, 0, 0);
            get_triplets_offset(N_cost, tripvec_P_i_nom, 0, Q_cost.cols());
            offset += Q_cost.rows();

            get_triplets_offset(R_cost, tripvec_P_i_nom, offset, offset);
            get_triplets_offset(N_cost.transpose(), tripvec_P_i_nom, offset, 0);
            offset += R_cost.rows();

            get_triplets_offset(Qu_ineq_cost, tripvec_P_i_nom, offset, offset);
            offset += Qu_ineq_cost.rows();

            P_i_nom.resize(offset, offset);
        }
        else
        {
            get_triplets_offset(Q_cost, tripvec_P_i_nom, 0, 0);
            get_triplets_offset(N_cost, tripvec_P_i_nom, 0, Q_cost.cols());
            offset += Q_cost.rows();

            get_triplets_offset(R_cost, tripvec_P_i_nom, offset, offset);
            get_triplets_offset(N_cost.transpose(), tripvec_P_i_nom, offset, 0);
            offset += R_cost.rows();

            P_i_nom.resize(offset, offset);
        }

        P_i_nom.setFromTriplets(tripvec_P_i_nom.begin(), tripvec_P_i_nom.end());
        q_i_nom = Eigen::VectorXd::Zero(P_i_nom.rows()); // init

        // reset off-nominal cost matrices
        P_i_vec.clear();
        q_i_vec.clear();

        // terminal cost
        // state vector: y = [x, sx_ineq]
        tripvec_P_N.reserve(P->nonZeros() + 
            std::max(Qx_term_ineq_cost.nonZeros(), Qx_ineq_cost.nonZeros()));
        
        if (terminal_state_ineq_specified)
        {
            if (soft_term_state_ineq)
            {
                get_triplets_offset(*P, tripvec_P_N, 0, 0);
                offset = P->rows();

                get_triplets_offset(Qx_term_ineq_cost, tripvec_P_N, offset, offset);
                offset += Qx_term_ineq_cost.rows();

                P_i_vec[n_horizon].resize(offset, offset);
            }
            else
            {
                get_triplets_offset(*P, tripvec_P_N, 0, 0);
                offset = P->rows();

                P_i_vec[n_horizon].resize(offset, offset);
            }
        }
        else if (soft_state_ineq)
        {
            get_triplets_offset(*P, tripvec_P_N, 0, 0);
            offset = P->rows();

            get_triplets_offset(Qx_ineq_cost, tripvec_P_N, offset, offset);
            offset += Qx_ineq_cost.rows();

            P_i_vec[n_horizon].resize(offset, offset);
        }
        else
        {
            get_triplets_offset(*P, tripvec_P_N, 0, 0);
            offset = P->rows();

            P_i_vec[n_horizon].resize(offset, offset);
        }

        P_i_vec[n_horizon].setFromTriplets(tripvec_P_N.begin(), tripvec_P_N.end());
        
        q_i_vec[n_horizon] = Eigen::VectorXd::Zero(P_i_vec[n_horizon].rows()); // init

    } // end if Hrep_constraints
    else
    {
        // nominal cost
        tripvec_P_i_nom.reserve(Q_cost.nonZeros() + R_cost.nonZeros() + 2*N_cost.nonZeros() + Zc_x.n + Zc_u.n);

        offset = 0;

        get_triplets_offset(Q_cost, tripvec_P_i_nom, 0, 0);
        get_triplets_offset(N_cost, tripvec_P_i_nom, 0, Q_cost.cols());
        offset += Q_cost.rows();

        get_triplets_offset(R_cost, tripvec_P_i_nom, offset, offset);
        get_triplets_offset(N_cost.transpose(), tripvec_P_i_nom, offset, 0);
        offset += R_cost.rows();

        offset += Zc_x.nG + Zc_u.nG;
        if (soft_state_ineq)
        {
            get_triplets_offset(Qx_ineq_cost, tripvec_P_i_nom, offset, offset);
            offset += Qx_ineq_cost.rows();
        }
        if (soft_input_ineq)
        {
            get_triplets_offset(Qu_ineq_cost, tripvec_P_i_nom, offset, offset);
            offset += Qu_ineq_cost.rows();
        }

        P_i_nom.resize(offset, offset);
        P_i_nom.setFromTriplets(tripvec_P_i_nom.begin(), tripvec_P_i_nom.end());
        q_i_nom = Eigen::VectorXd::Zero(P_i_nom.rows()); // init

        // reset off-nominal cost matrices
        P_i_vec.clear();
        q_i_vec.clear();

        // terminal cost
        tripvec_P_N.reserve(P->nonZeros() + 
            std::max(Qx_term_ineq_cost.nonZeros(), Qx_ineq_cost.nonZeros()));

        offset = 0;
        get_triplets_offset(*P, tripvec_P_N, offset, offset);
        offset += P->rows();

        if (terminal_state_ineq_specified)
            offset += Zc_x_term.nG;
        else
            offset += Zc_x.nG;

        Eigen::SparseMatrix<double> * Qx_term_ptr = nullptr;
        Eigen::SparseMatrix<double> empty_mat (0,0);
        if (terminal_state_ineq_specified && soft_term_state_ineq)
            Qx_term_ptr = &Qx_term_ineq_cost;
        else if (soft_state_ineq)
            Qx_term_ptr = &Qx_ineq_cost;
        else
            Qx_term_ptr = &empty_mat;

        get_triplets_offset(*Qx_term_ptr, tripvec_P_N, offset, offset);
        offset += Qx_term_ptr->rows();

        P_i_vec[n_horizon].resize(offset, offset);
        P_i_vec[n_horizon].setFromTriplets(tripvec_P_N.begin(), tripvec_P_N.end());
        q_i_vec[n_horizon] = Eigen::VectorXd::Zero(P_i_vec[n_horizon].rows()); // init
    }

    // make gradient
    make_cost_gradient_from_reference();
}

// reference tracking
void MPC_QP::make_cost_gradient_from_reference()
{
    // stage cost
    auto it_x_ref = x_ref.begin();
    for (int k=1; k<=n_horizon; ++k)
    {
        // create q_i_vec[k] if doesn't already exist
        if (!q_i_vec.count(k))
            q_i_vec[k] = q_i_nom;

        // reference contribution to gradient
        if (k == n_horizon && terminal_cost_specified)
        {
            q_i_vec[k].segment(0, n_state) = q_xN_cost - P_cost * (*it_x_ref);
        }
        else
        {
            q_i_vec[k].segment(0, n_state) = q_x_cost - Q_cost * (*it_x_ref);
            q_i_vec[k].segment(n_state, n_input) = q_u_cost;
        }

        // increment reference
        ++it_x_ref;
    }
}

void MPC_QP::make_const_cost_from_reference()
{
    // init const cost
    b = 0;

    // loop through to add reference tracking cost
    for (int k=1; k<n_horizon; ++k)
    {
        // add reference tracking cost
        b += x_ref[k-1].transpose() * Q_cost * x_ref[k-1];
    }
    b += x_ref.back().transpose() * P_cost * x_ref.back();
    b /= 2;
}

// check problem validity
std::pair<bool, std::string> MPC_QP::check_problem_validity()
{
    // declare
    Eigen::SparseMatrix<double> *P_k = nullptr, *C_k = nullptr, *D_k = nullptr, *D_kp1 = nullptr, *G_k = nullptr;
    Eigen::VectorXd *q_k = nullptr, *crhs_k = nullptr, *w_k = nullptr; 
    std::stringstream ss;
    int n_y, n_eq, n_ineq;
    bool valid = true;

    // loop through time steps
    for (int k=0; k<=n_horizon; k++)
    {
        // get matrices / vectors
        P_k = (P_i_vec.count(k)) ? &P_i_vec[k] : &P_i_nom;
        q_k = (q_i_vec.count(k)) ? &q_i_vec[k] : &q_i_nom;
        C_k = (C_i_vec.count(k)) ? &C_i_vec[k] : &C_i_nom;
        D_k = (D_i_vec.count(k)) ? &D_i_vec[k] : &D_i_nom;
        D_kp1 = (D_i_vec.count(k+1)) ? &D_i_vec[k+1] : &D_i_nom;
        crhs_k = (crhs_i_vec.count(k)) ? &crhs_i_vec[k] : &crhs_i_nom;
        G_k = (G_i_vec.count(k)) ? &G_i_vec[k] : &G_i_nom;
        w_k = (w_i_vec.count(k)) ? &w_i_vec[k] : &w_i_nom;

        // check state dimensions
        n_y = n_y_vec[k];

        if (P_k->rows() != n_y || P_k->cols() != n_y)
        {
            ss << "Error: P matrix at time step " << k << " has inconsistent dimensions." << std::endl;
            valid = false;
        }
        if (q_k->size() != n_y)
        {
            ss << "Error: q vector at time step " << k << " has inconsistent dimensions." << std::endl;
            valid = false;
        }
        if (k != n_horizon && C_k->cols() != n_y)
        {
            ss << "Error: C matrix at time step " << k << " has inconsistent state dimensions." << std::endl;
            valid = false;
        }
        if (k != 0 && D_k->cols() != n_y)
        {
            ss << "Error: D matrix at time step " << k << " has inconsistent state dimensions." << std::endl;
            valid = false;
        }
        if (G_k->cols() != n_y)
        {
            ss << "Error: G matrix at time step " << k << " has inconsistent state dimensions." << std::endl;
            valid = false;
        }
        
        // check equality constraint dimensions
        if (k != n_horizon)
        {
            n_eq = n_eq_vec[k];
            if (C_k->rows() != n_eq)
            {
                ss << "Error: C matrix at time step " << k << " has inconsistent equality constraint dimensions." << std::endl;
                valid = false;
            }
            if (D_kp1->rows() != n_eq)
            {
                ss << "Error: D matrix at time step " << k+1 << " has inconsistent equality constraint dimensions." << std::endl;
                valid = false;
            }
            if (crhs_k->size() != n_eq)
            {
                ss << "Error: crhs vector at time step " << k << " has inconsistent equality constraint dimensions." << std::endl;
                valid = false;
            }
        }

        // check inequality constraint dimensions
        n_ineq = n_ineq_vec[k];

        if (G_k->rows() != n_ineq)
        {
            ss << "Error: G matrix at time step " << k << " has inconsistent inequality constraint dimensions." << std::endl;
            valid = false;
        }
        if (w_k->size() != n_ineq)
        {
            ss << "Error: w vector at time step " << k << " has inconsistent inequality constraint dimensions." << std::endl;
            valid = false;
        }
    }

    return std::make_pair(valid, ss.str());
}

// get triplets for matrix
void MPC_QP::get_triplets_offset(const Eigen::SparseMatrix<double> &mat, std::vector<Eigen::Triplet<double>> &triplets, 
            int i_offset, int j_offset)
{
    // get triplets
    for (int k=0; k<mat.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
            triplets.push_back(Eigen::Triplet<double>(it.row() + i_offset, it.col() + j_offset, it.value()));
        }
    }
}

void MPC_QP::get_triplets_offset(const Eigen::Ref<const Eigen::VectorXd> vec, std::vector<Eigen::Triplet<double>> &triplets, 
            int i_offset, int j_offset)
{
    // get triplets
    for (int i=0; i<vec.size(); ++i)
    {
        triplets.push_back(Eigen::Triplet<double>(i + i_offset, j_offset, vec(i)));
    }
}

void MPC_QP::get_triplets_offset_identity(int n_dim, std::vector<Eigen::Triplet<double>> &triplets, int i_offset, int j_offset)
{
    // get triplets
    for (int i=0; i<n_dim; ++i)
    {
        triplets.push_back(Eigen::Triplet<double>(i + i_offset, i + j_offset, 1.0));
    }
}

void MPC_QP::get_triplets_offset_neg_identity(int n_dim, std::vector<Eigen::Triplet<double>> &triplets, int i_offset, int j_offset)
{
    // get triplets
    for (int i=0; i<n_dim; ++i)
    {
        triplets.push_back(Eigen::Triplet<double>(i + i_offset, i + j_offset, -1.0));
    }
}

// get indices
std::vector<int> MPC_QP::get_indices(int offset, int range)
{
    std::vector<int> indices(range);
    std::iota(indices.begin(), indices.end(), offset);
    return indices;
}