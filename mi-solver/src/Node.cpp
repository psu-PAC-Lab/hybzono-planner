#include "Node.hpp"

using namespace MI_QPIPMPC;

// constructor
Node::Node(const Data& data, const QP_IP_MPC::QP_primal_dual_mpc &solver,
           double lower, const std::map<int, std::vector<int>> &region_vec)
    : solver(solver), lower(lower), region_vec(region_vec)
{
    // get depth
    get_depth(data);

    // init other member variables
    frac_idx = 0;
    num_iter = 0;
    qp_solve_time = 0;
    status = QP_NO_SOL;
    frac = 0;
    cons_violated = true;
    
    // node is valid
    valid = true;
}

Node::Node(const Data& data, const QP_IP_MPC::QP_primal_dual_mpc &solver,
           double lower, const std::map<int, std::vector<int>> &region_vec,
           const std::vector<int> &term_region_vec)
    : solver(solver), lower(lower), region_vec(region_vec), term_region_vec(term_region_vec)
{
    // get depth
    get_depth(data);
    get_terminal_constraint_depth();

    // terminal constraints enabled
    terminal_constraint = true;

    // init other member variables
    frac_idx = 0;
    num_iter = 0;
    qp_solve_time = 0;
    status = QP_NO_SOL;
    frac = 0;
    cons_violated = true;
    
    // node is valid
    valid = true;
}

// sorting methods - returns true if this worse that n
bool Node::operator < (const Node &n) const
{
    if (n.priority)
        return true;
    else if (this->priority)
        return false;
    else if (this->lower != n.lower)
        return this->lower > n.lower;
    else if (this->frac != n.frac)
        return this->frac < n.frac;
    else
        return this->term_frac < n.term_frac;
}

double Node::operator - (const Node &n) const
{
    return this->lower - n.lower;
}

// set methods
void Node::set_region_vec(const std::map<int, std::vector<int>> &region_vec)
{
    this->region_vec = region_vec;
}

void Node::set_frac(double frac)
{
    this->frac = frac;
}

void Node::set_term_frac(double term_frac)
{
    this->term_frac = term_frac;
}

// get methods
void Node::get_depth(const Data& data)
{
    depth = 0;
    for (int k=0; k<=data.n_horizon; k++)
    {
        if (region_vec.count(k) && (region_vec[k].size() == 1))
            depth++;
    }
}

void Node::get_terminal_constraint_depth()
{
    if (term_region_vec.size() <= 1)
        term_constraint_depth = 1;
    else
        term_constraint_depth = 0;
}

// reset solver
void Node::reset_solver(const Data& data)
{
    // re-initialize solver
    solver.setup(data.qpipmpc_data.P_i_nom, data.qpipmpc_data.q_i_nom,
        data.qpipmpc_data.C_i_nom, data.qpipmpc_data.D_i_nom,
        data.qpipmpc_data.crhs_i_nom, data.qpipmpc_data.G_i_nom,
        data.qpipmpc_data.w_i_nom, data.n_horizon, data.qpipmpc_data.b);

    solver.set_P_vec(data.qpipmpc_data.P_i_vec);
    solver.set_q_vec(data.qpipmpc_data.q_i_vec);
    solver.set_C_vec(data.qpipmpc_data.C_i_vec);
    solver.set_D_vec(data.qpipmpc_data.D_i_vec);
    solver.set_crhs_vec(data.qpipmpc_data.crhs_i_vec);
    solver.set_G_vec(data.qpipmpc_data.G_i_vec);
    solver.set_w_vec(data.qpipmpc_data.w_i_vec);
    solver.set_settings(data.qp_settings);
}

// solve
void Node::solve(const Data& data)
{
    // update constraints
    solver.set_P_vec(P_i_vec);
    solver.set_q_vec(q_i_vec);
    solver.set_C_vec(C_i_vec);
    solver.set_D_vec(D_i_vec);
    solver.set_crhs_vec(crhs_i_vec);
    solver.set_G_vec(G_i_vec);
    solver.set_w_vec(w_i_vec);

    // solve problem
    qp_results = solver.solve();
    qp_solve_time = qp_results.sol_time;
    x = qp_results.y;
    num_iter = qp_results.num_iter;

    // map solution back to original problem
    add_zeroed_variables(data, x);

    // set solver status
    if (!qp_results.feas || qp_results.objective > QP::inf)
        status = QP_INFEASIBLE;
    else if (qp_results.converged)
        status = QP_SOLVED;
    else
        status = QP_SUBOPTIMAL;

    // update objective function lower bound if solver is converged
    switch (status)
    {
        case QP_SOLVED: case QP_SUBOPTIMAL:
        {
            obj = qp_results.objective;
            lower = qp_results.objective;
            break;
        }
        case QP_INFEASIBLE:
        {
            obj = QP::inf;
            lower = QP::inf;
            break;
        }
    }
}

// check for H-rep constraint violations at given time step
std::pair<bool, int> Node::check_Hrep_constraint_violation(int k_timestep, const Data& data)
{
    // get state solution and position at step k
    Eigen::VectorXd pos_k = data.qpipmpc_data.Pc*x(data.qpipmpc_data.idx_state.at(k_timestep));

    // check whether there is a region without constraint violation and get its index
    double max_err = QP::inf; // init
    int region_no_violation = -1; // init
    bool constraint_violation = true; // init
    double err; // declare
    for (const int &region : region_vec.at(k_timestep))
    {
        err = (data.A_Hrep.at(region)*pos_k - data.b_Hrep.at(region)).maxCoeff(); // invalid matrix product
        if (err < max_err)
        {
            max_err = err;
            region_no_violation = region;
        }
    }

    if (max_err <= data.settings.eps_feas)
        constraint_violation = false;

    // return constraint violation, region_no_violation
    return std::pair<bool, int> {constraint_violation, region_no_violation};
}

// get regions without constraint violation
void Node::get_regions_no_constraint_violation(const Data& data)
{
    // init
    region_vec_no_violation.clear();
    cons_violated = false;

    // loop through by time step and check for constraint violations
    for (int k=0; k<=data.n_horizon; k++)
    {
        // handle case where region already specified
        if (region_vec.count(k) && (region_vec[k].size() == 1))
            region_vec_no_violation[k] = { (*(region_vec[k].begin())) };
        else
        {
            // check for constraint violation and get index of
            // region w/o constraint violation if it exists
            std::pair<bool, int> cons_check = check_Hrep_constraint_violation(k, data);
            if (!cons_check.first)
                region_vec_no_violation[k] = { cons_check.second };
            else // constraints violated
            {
                k_first_cons_violation = k;
                cons_violated = true;
                break;
            }
        }
    }
}

// get fractionality vector
Eigen::VectorXd Node::get_frac_vec(int k, const std::vector<int> &regions, const Data& data) const
{
    // get indices of binary (region) variables within node solution
    // at time step k and corresponding to input region vector
    std::vector<int> idx_regions;
    for (auto it = regions.begin(); it != regions.end(); ++it)
        idx_regions.push_back(data.qpipmpc_data.idx_binvar.at(k)[*it-1]); // -1 to correct for one based indexing

    // return fractional values of binary variables
    return x(idx_regions);
}

// terminal constraint fractionality vector
Eigen::VectorXd Node::get_term_frac_vec(const std::vector<int> &regions, const Data& data) const
{
    // make sure terminal constraint is enabled
    if (!this->terminal_constraint) return Eigen::VectorXd::Zero(0);

    // get indices of binary (region) variables within node solution
    // at time step k and corresponding to input region vector
    std::vector<int> idx_regions;
    for (auto it = regions.begin(); it != regions.end(); ++it)
        idx_regions.push_back(data.idx_binvar_term[*it-1]); // -1 to correct for one based indexing

    // return fractional values of binary variables
    return x(idx_regions);
}

// reachability logic
bool Node::init_reachability_constraints(const Data& data) 
{
    // time step 0
    region_vec[0].push_back(data.region_0);

    // x0 to region reachability
    for (int k=1; k<=data.n_horizon; k++)
        region_vec[k] = data.R_x02region.at(k);

    // make consistent
    make_reachability_constraints_consistent(data, region_vec);

    // check valid
    bool valid = check_reachability_valid(data);

    // apply constraints
    if (valid)
        apply_reachability_constraints(data);

    return valid;
}

bool Node::update_reachability_constraints(const Data &data, int kb, 
    const std::vector<int> &regions_kb) 
{
    region_vec[kb] = regions_kb;

    // make consistent
    make_reachability_constraints_consistent(data, region_vec);

    // check valid
    bool valid = check_reachability_valid(data);
    
    // apply constraints
    if (valid)
        apply_reachability_constraints(data);
    
    return valid;
}


bool Node::check_reachability_valid(const Data& data)
{
    // check reachability constraints aren't empty at any time points
    for (int k=0; k<=data.n_horizon; k++)
    {
        if (region_vec[k].empty())
            return false;
    }

    return true;
}


void Node::update_terminal_constraint(const std::vector<int> &term_region_vec)
{
    this->term_region_vec = term_region_vec;
    terminal_constraint = true;
}

void Node::make_reachability_constraints_consistent(const Data& data, std::map<int, std::vector<int>>& region_vec)
{
    // working variable: region to region reachability
    std::vector<int> regions_r2r;
    regions_r2r.reserve(data.n_regions*(data.n_horizon+1));

    // boolean array declarations and lambdas
    bool * regions_k_arr = new bool [data.n_regions];
    bool * regions_r2r_arr = new bool [data.n_regions];

    auto set_false = [&] (bool * arr, int len)
    {
        for (int i=0; i<len; i++)
            arr[i] = false;
    };

    auto bitwise_and = [&] (bool * arr1, bool * arr2, int len)
    {
        for (int i=0; i<len; i++)
            arr1[i] = arr1[i] && arr2[i];
    };

    auto bitwise_or_vec = [&] (bool * arr1, const std::vector<int> &vec2)
    {
        for (auto it = vec2.begin(); it != vec2.end(); ++it)
            arr1[*it-1] = true; // zero-based indexing
    };

    auto set_from_vector = [&] (bool * arr, const std::vector<int> &vec)
    {
        set_false(arr, data.n_regions);
        for (auto it = vec.begin(); it != vec.end(); ++it)
            arr[*it-1] = true; // zero-based indexing
    };

    auto set_vector_from_arr = [&] (std::vector<int> &vec, bool * arr, int len)
    {
        vec.clear();
        for (int i=0; i<len; i++)
        {
            if (arr[i])
                vec.push_back(i+1); // one-based indexing
        }
    };

    // loop through by time step
    int k, kr; // declare
    for (auto it_leaf_k = region_vec.begin(); it_leaf_k != region_vec.end(); ++it_leaf_k)
    {
        // get time step
        k = it_leaf_k->first;

        // no computation necessary if region_vec[k] has only one element
        if (it_leaf_k->second.size() == 1) 
            continue;          

        // get binary variable array corresponding for valid regions at time step k
        set_from_vector(regions_k_arr, it_leaf_k->second);

        // region to region reachability
        for (auto it_leaf_kr = region_vec.begin(); it_leaf_kr != region_vec.end(); ++it_leaf_kr)
        {
            // get time step
            kr = it_leaf_kr->first;

            // reset regions_r2r
            set_false(regions_r2r_arr, data.n_regions);

            // skip if kr = k since this would be zero step reachability
            if (kr == k) continue;

            // get possible regions at time step k given possible regions at time step kr
            const std::map<int, std::vector<int>> * R_region2region_ptr = &(data.R_region2region.at(abs(k-kr))); // only do this map lookup once
            for (auto it_kr = it_leaf_kr->second.begin(); it_kr != it_leaf_kr->second.end(); ++it_kr)
            {
                bitwise_or_vec(regions_r2r_arr, R_region2region_ptr->at(*it_kr));
            }

            // apply region-to-region reachability to regions at time step k
            bitwise_and(regions_k_arr, regions_r2r_arr, data.n_regions);
        }

        // apply reachability constraints
        set_vector_from_arr(it_leaf_k->second, regions_k_arr, data.n_regions);
    }

    // free memory
    delete [] regions_k_arr;
    delete [] regions_r2r_arr;
}

void Node::apply_reachability_constraints(const Data& data)
{
    // init / declarations
    int idx_node = 0, idx_prob = 0;
    idx_node2problem.clear();
    std::vector<int> idx_zeroed;
    idx_zeroed.reserve(data.n_regions);
    std::vector<int> zero_rows_common;
    const Eigen::SparseMatrix<double> * C_ptr = nullptr, * D_ptr = nullptr;
    const Eigen::VectorXd * crhs_ptr = nullptr;
    std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> CD_k_ptr;
    std::vector<int> CD_k_cols;
    CD_k_cols.reserve(2*data.qpipmpc_data.idx_y.at(0).size());

    // loop through time steps
    for (int k=0; k <= data.n_horizon; k++)
    {
        // get indices of binary variables
        const std::vector<int> * idx_ptr = nullptr;
        if (data.qpipmpc_data.idx_i_vec.count(k))
            idx_ptr = &(data.qpipmpc_data.idx_i_vec.at(k));
        else
            idx_ptr = &(data.qpipmpc_data.idx_i_nom);

        // get indices of zeroed binary variables
        idx_zeroed.clear();
        auto it_region_vec = region_vec[k].begin(); // must be sorted
        for (int region = 1; region <= data.n_regions; region++)
        {
            if (it_region_vec != region_vec[k].end() && region == *it_region_vec)
                ++it_region_vec;
            else
                idx_zeroed.push_back((*idx_ptr)[region-1]); // zero-based indexing
        }

        // check for continuous variables that are forced to zero by zeroed binaries
        // and add to idx_zeroed
        if (k < data.n_horizon)
        {
            C_ptr = data.qpipmpc_data.C_i_vec.count(k) ? &data.qpipmpc_data.C_i_vec.at(k) : &data.qpipmpc_data.C_i_nom;
            D_ptr = data.qpipmpc_data.D_i_vec.count(k+1) ? &data.qpipmpc_data.D_i_vec.at(k+1) : &data.qpipmpc_data.D_i_nom;
            crhs_ptr = data.qpipmpc_data.crhs_i_vec.count(k) ? &data.qpipmpc_data.crhs_i_vec.at(k) : &data.qpipmpc_data.crhs_i_nom;
            hor_cat_row_major(C_ptr, D_ptr, CD_k_ptr);
        }
        // else can reuse CD_k_ptr from k-1

        for (int i=0; i<CD_k_ptr->outerSize(); i++)
        {
            bool row_zeroed = true; // init
            CD_k_cols.clear();
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(*CD_k_ptr, i); it; ++it)
            {
                if (k < data.n_horizon && (*crhs_ptr)[it.row()] == 0 && it.col() < C_ptr->cols())
                {
                    if ((it.value() > 0 &&
                        std::find(data.idx_pos_vars.at(k).begin(), data.idx_pos_vars.at(k).end(), it.col()) != data.idx_pos_vars.at(k).end()) ||
                        std::find(idx_zeroed.begin(), idx_zeroed.end(), it.col()) != idx_zeroed.end())
                    {
                        CD_k_cols.push_back(it.col());
                    }
                    else
                    {
                        row_zeroed = false;
                    }
                }
                else if (k == data.n_horizon && (*crhs_ptr)[it.row()] == 0 && it.col() >= C_ptr->cols())
                {
                    if ((it.value() > 0 && 
                        std::find(data.idx_pos_vars.at(k).begin(), data.idx_pos_vars.at(k).end(), it.col()-C_ptr->cols()) != data.idx_pos_vars.at(k).end()) ||
                        std::find(idx_zeroed.begin(), idx_zeroed.end(), it.col()-C_ptr->cols()) != idx_zeroed.end())
                    {
                        CD_k_cols.push_back(it.col()-C_ptr->cols());
                    }
                    else
                    {
                        row_zeroed = false;
                    }
                }
                else
                {
                    row_zeroed = false;
                }
                if (!row_zeroed) break;
            }
            if (row_zeroed)
            {
                for (int j=0; j<CD_k_cols.size(); j++)
                {
                    if (std::find(idx_zeroed.begin(), idx_zeroed.end(), CD_k_cols[j]) == idx_zeroed.end())
                        idx_zeroed.push_back(CD_k_cols[j]);
                }
            }
        }
        std::sort(idx_zeroed.begin(), idx_zeroed.end());

        // track mapping of variables in sub-problem to original problem
        auto it_idx_zeroed = idx_zeroed.begin();
        for (int idx_prob_k = 0; idx_prob_k < data.qpipmpc_data.idx_y.at(k).size(); idx_prob_k++)
        {
            if (it_idx_zeroed != idx_zeroed.end() && idx_prob_k == *it_idx_zeroed)
            {
                ++it_idx_zeroed;
            }
            else
            {
                idx_node2problem.push_back(std::make_pair(idx_node, idx_prob));
                idx_node++;
            }

            idx_prob++;
        }

        // remove zeroed variables from cost
        P_i_vec[k] = data.qpipmpc_data.P_i_vec.count(k) ? data.qpipmpc_data.P_i_vec.at(k) : data.qpipmpc_data.P_i_nom;
        q_i_vec[k] = data.qpipmpc_data.q_i_vec.count(k) ? data.qpipmpc_data.q_i_vec.at(k) : data.qpipmpc_data.q_i_nom;

        remove_rows_cols(&P_i_vec[k], idx_zeroed, idx_zeroed);
        remove_elements(&q_i_vec[k], idx_zeroed);

        // remove zeroed binary variables from equality constraints
        if (k != data.n_horizon)
        {
            C_i_vec[k] = data.qpipmpc_data.C_i_vec.count(k) ? data.qpipmpc_data.C_i_vec.at(k) : data.qpipmpc_data.C_i_nom;
            crhs_i_vec[k] = data.qpipmpc_data.crhs_i_vec.count(k) ? data.qpipmpc_data.crhs_i_vec.at(k) : data.qpipmpc_data.crhs_i_nom;

            remove_cols(&C_i_vec[k], idx_zeroed);
        }
        if (k != 0)
        {
            D_i_vec[k] = data.qpipmpc_data.D_i_vec.count(k) ? data.qpipmpc_data.D_i_vec.at(k) : data.qpipmpc_data.D_i_nom;

            remove_cols(&D_i_vec[k], idx_zeroed);
        }

        // remove zeroed binary variables from inequality constraints
        G_i_vec[k] = data.qpipmpc_data.G_i_vec.count(k) ? data.qpipmpc_data.G_i_vec.at(k) : data.qpipmpc_data.G_i_nom;
        w_i_vec[k] = data.qpipmpc_data.w_i_vec.count(k) ? data.qpipmpc_data.w_i_vec.at(k) : data.qpipmpc_data.w_i_nom;

        remove_cols(&G_i_vec[k], idx_zeroed);
    }

    // check for constraints that may have been removed as a result of deleting variables
    for (int k=0; k<=data.n_horizon; k++)
    {
        // equality constraints
        if (k != data.n_horizon)
        {
            std::vector<int> zero_rows_C = get_zero_rows(&C_i_vec[k]);
            std::vector<int> zero_rows_D = get_zero_rows(&D_i_vec[k+1]);

            zero_rows_common.clear();
            for (auto it = zero_rows_C.begin(); it != zero_rows_C.end(); ++it)
            {
                if (std::find(zero_rows_D.begin(), zero_rows_D.end(), *it) != zero_rows_D.end())
                    zero_rows_common.push_back(*it);
            }

            if (!zero_rows_common.empty())
            {
                remove_rows(&C_i_vec[k], zero_rows_common);
                remove_rows(&D_i_vec[k+1], zero_rows_common);
                remove_elements(&crhs_i_vec[k], zero_rows_common);
            }
        }

        // inequality constraints
        std::vector<int> zero_rows_G = get_zero_rows(&G_i_vec[k]);

        if (!zero_rows_G.empty())
        {
            remove_rows(&G_i_vec[k], zero_rows_G);
            remove_elements(&w_i_vec[k], zero_rows_G);
        }
    }

    // apply terminal constraints if applicable
    if (this->terminal_constraint)
        apply_terminal_constraint(data);
}

// add zeroed variables back to solution
void Node::add_zeroed_variables(const Data& data, Eigen::VectorXd &x)
{
    // total size of x
    int n_y_tot = 0;
    for (auto it = data.qpipmpc_data.idx_y.begin(); it != data.qpipmpc_data.idx_y.end(); ++it)
        n_y_tot += it->second.size();

    // set to zero
    Eigen::VectorXd x_new (n_y_tot);
    x_new.setZero();

    // map node solution to original problem
    for (auto it = idx_node2problem.begin(); it != idx_node2problem.end(); ++it)
        x_new(it->second) = x(it->first);

    // output
    x = x_new;
}


// terminal constraint logic
void Node::apply_terminal_constraint(const Data& data)
{
    // make sure terminal constraint is enabled
    if (!this->terminal_constraint) return;

    // conditional to prevent duplicate rows
    if (term_region_vec.size() == data.n_term_regions) return;

    // matrix dimensions
    int m_C = C_i_vec[data.n_horizon-1].rows();
    int n_C = C_i_vec[data.n_horizon-1].cols();
    int m_D = D_i_vec[data.n_horizon].rows();
    int n_D = D_i_vec[data.n_horizon].cols();
    int m_crhs = crhs_i_vec[data.n_horizon-1].size();

    // resize matrices
    C_i_vec[data.n_horizon-1].conservativeResize(m_C+1, n_C);
    D_i_vec[data.n_horizon].conservativeResize(m_D+1, n_D);
    crhs_i_vec[data.n_horizon-1].conservativeResize(m_crhs+1);

    // get indices of terminal binary variables for reduced-variable node
    std::vector<int> idx_term_node;
    idx_term_node.reserve(data.idx_term.size());
    for (auto it = data.idx_term.begin(); it != data.idx_term.end(); ++it)
    {
        for (auto it_node2problem = idx_node2problem.begin(); it_node2problem != idx_node2problem.end(); ++it_node2problem)
        {
            if (it_node2problem->second == *it)
            {
                idx_term_node.push_back(it_node2problem->first);
                break;
            }
        }
    }

    if (idx_term_node.size() != data.n_term_regions)
        throw std::runtime_error("Terminal constraint size mismatch");

    // add terminal constraint selection
    for (auto it = term_region_vec.begin(); it != term_region_vec.end(); ++it)
        D_i_vec[data.n_horizon].insert(m_D, idx_term_node[*it-1]) = 1; // zero-based indexing
    crhs_i_vec[data.n_horizon-1](m_crhs) = -1; // zero-based indexing
}

// give priority to node
void Node::give_priority()
{
    priority = true;
}

void Node::remove_priority()
{
    priority = false;
}

bool Node::is_priority() const
{
    return priority;
}

// check if valid
bool Node::is_valid() const
{
    return valid;
}

// get approximate cost of trajectory + region_vec
double Node::get_approx_cost(const Data& data, const std::map<int, std::vector<int>> &region_vec) const
{
    // init
    double cost = 0;
    Eigen::VectorXd x_app = this->x;    
    Eigen::VectorXd x_bin_k, x_app_k; // declare
    const Eigen::SparseMatrix<double> * P_k = nullptr;
    const Eigen::VectorXd * q_k = nullptr;

    // loop through time steps
    for (int k=0; k<=data.n_horizon; k++)
    {
        // should only be one region per time step
        if (region_vec.at(k).size() > 1) 
        {
            cost = QP::inf;
            return cost;
        }
        else
        {
            // binary variables
            x_bin_k = Eigen::VectorXd::Zero(data.qpipmpc_data.idx_binvar.at(k).size());
            x_bin_k(region_vec.at(k)[0]-1) = 1; // 1-based indexing -> 0-based indexing

            // insert into solution
            x_app(data.qpipmpc_data.idx_binvar.at(k)) = x_bin_k;
        }

        // get cost of time step k
        if (data.qpipmpc_data.P_i_vec.count(k))
            P_k = &(data.qpipmpc_data.P_i_vec.at(k));
        else
            P_k = &(data.qpipmpc_data.P_i_nom);

        if (data.qpipmpc_data.q_i_vec.count(k))
            q_k = &(data.qpipmpc_data.q_i_vec.at(k));
        else
            q_k = &(data.qpipmpc_data.q_i_nom);

        x_app_k = x_app(data.qpipmpc_data.idx_y.at(k));
        cost += 0.5*(x_app_k.transpose()).dot((*P_k)*x_app_k) + (*q_k).dot(x_app_k);            
    }

    // return cost
    return cost;
}

// utilities
void Node::remove_cols(Eigen::SparseMatrix<double> * M, const std::vector<int> &cols)
{
    if (!std::is_sorted(cols.begin(), cols.end()))
        throw std::invalid_argument("cols must be sorted");

    // check for early exit
    if (cols.empty())
        return;

    M->makeCompressed();

    // init
    std::vector<Eigen::Triplet<double>> tripvec;
    tripvec.reserve(M->nonZeros());
    
    int col;
    int ind_cols = 0;
    int col_active = -1;
    int col_next = cols[ind_cols];

    // loop through M
    for (int k=0; k<M->outerSize(); ++k)
    {
        col = k;
        if (col == col_next)
        {
            // update column
            col_active = col_next;
            ind_cols++;
            if (ind_cols < cols.size())
                col_next = cols[ind_cols];
            else
                col_next = -1;
            
            // skip
            continue;
        }

        for (Eigen::SparseMatrix<double>::InnerIterator it(*M, k); it; ++it)
        {           
            // add triplet
            tripvec.push_back(Eigen::Triplet<double>(it.row(), col-ind_cols, it.value()));
        }
    }

    // update matrix
    M->resize(M->rows(), M->cols()-cols.size());
    M->setFromSortedTriplets(tripvec.begin(), tripvec.end());
}

void Node::remove_rows_cols(Eigen::SparseMatrix<double> * M, const std::vector<int> &rows, const std::vector<int> &cols)
{
    if (!std::is_sorted(cols.begin(), cols.end()))
        throw std::invalid_argument("cols must be sorted");
    if (!std::is_sorted(rows.begin(), rows.end()))
        throw std::invalid_argument("rows must be sorted");

    M->makeCompressed();

    // check for early exit
    if (cols.empty() && rows.empty())
        return;

    // init
    std::vector<Eigen::Triplet<double>> tripvec;
    tripvec.reserve(M->nonZeros());
    
    int col, row;
    int ind_cols = 0, ind_rows;
    int col_next = cols[ind_cols];

    // loop through M
    for (int k=0; k<M->outerSize(); ++k)
    {
        col = k;
        if (col == col_next)
        {
            // update column
            ind_cols++;
            if (ind_cols < cols.size())
                col_next = cols[ind_cols];
            else
                col_next = -1;
            
            // skip
            continue;
        }

        for (Eigen::SparseMatrix<double>::InnerIterator it(*M, k); it; ++it)
        {
            // find new row
            row = it.row();
            ind_rows = 0;
            while (ind_rows < rows.size() && row > rows[ind_rows])
                ind_rows++;

            if (ind_rows < rows.size() && row == rows[ind_rows])
                continue;

            // add triplet
            tripvec.push_back(Eigen::Triplet<double>(row-ind_rows, col-ind_cols, it.value()));
        }
    }

    // update matrix
    M->resize(M->rows()-rows.size(), M->cols()-cols.size());
    M->setFromSortedTriplets(tripvec.begin(), tripvec.end());
}

void Node::remove_rows(Eigen::SparseMatrix<double> * M, const std::vector<int> &rows)
{
    if (!std::is_sorted(rows.begin(), rows.end()))
        throw std::invalid_argument("rows must be sorted");

    // check for early exit
    if (rows.empty())
        return;

    // init
    std::vector<Eigen::Triplet<double>> tripvec;
    tripvec.reserve(M->nonZeros());
    
    int row, ind_rows;

    // loop through M
    for (int k=0; k<M->outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(*M, k); it; ++it)
        {            
            // find new row
            row = it.row();
            ind_rows = 0;
            while (ind_rows < rows.size() && row > rows[ind_rows])
                ind_rows++;

            if (ind_rows < rows.size() && row == rows[ind_rows])
                continue;
            
            // add triplet
            tripvec.push_back(Eigen::Triplet<double>(row-ind_rows, it.col(), it.value()));
        }
    }

    // update matrix
    M->resize(M->rows()-rows.size(), M->cols());
    M->setFromSortedTriplets(tripvec.begin(), tripvec.end());
}

void Node::remove_elements(Eigen::VectorXd * v, const std::vector<int> &elems)
{
    if (!std::is_sorted(elems.begin(), elems.end()))
        throw std::invalid_argument("rows must be sorted");

    // check for early exit
    if (elems.empty())
        return;

    // init
    Eigen::VectorXd v_new (v->size()-elems.size());

    // populate
    auto it_elems = elems.begin();
    int i_new = 0;
    for (int i=0; i<v->size(); i++)
    {
        if (it_elems != elems.end() && i == *it_elems)
        {
            ++it_elems;
            continue;
        }    
        else
        {
            v_new(i_new) = (*v)(i);
            i_new++;
        }
    }

    // update
    *v = v_new;
}

std::vector<int> Node::get_zero_rows(const Eigen::SparseMatrix<double> * M)
{
    // declare
    std::vector<int> zero_rows, nonzero_rows;
    zero_rows.reserve(M->rows());
    nonzero_rows.reserve(M->rows());

    // get nonzero rows
    for (int k=0; k<M->outerSize(); k++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(*M, k); it; ++it)
        {
            if (nonzero_rows.size() == M->rows()) // early exit
                return zero_rows;

            if (std::find(nonzero_rows.begin(), nonzero_rows.end(), it.row()) == nonzero_rows.end())
                nonzero_rows.push_back(it.row());
        }
    }
    std::sort(nonzero_rows.begin(), nonzero_rows.end());

    // get zero rows
    auto it = nonzero_rows.begin();
    for (int i=0; i<M->rows(); i++)
    {
        if (it != nonzero_rows.end() && i == *it)
            ++it;
        else
            zero_rows.push_back(i);
    }

    return zero_rows;
}

// get indices
std::vector<int> Node::get_indices(int offset, int range)
{
    std::vector<int> indices(range);
    std::iota(indices.begin(), indices.end(), offset);
    return indices;
}

// horizontal contatenation into a row major matrix
void Node::hor_cat_row_major(const Eigen::SparseMatrix<double> * A_ptr, const Eigen::SparseMatrix<double> * B_ptr, std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> &AB_ptr)
{
    Eigen::SparseMatrix<double, Eigen::RowMajor> A_rm = *A_ptr;
    Eigen::SparseMatrix<double, Eigen::RowMajor> B_rm = *B_ptr;

    AB_ptr.reset();
    AB_ptr = std::make_shared<Eigen::SparseMatrix<double, Eigen::RowMajor>>(A_ptr->rows(), A_ptr->cols() + B_ptr->cols());
    AB_ptr->reserve(A_rm.nonZeros() + B_rm.nonZeros());

    for (Eigen::Index k=0; k<A_rm.outerSize(); ++k)
    {
        AB_ptr->startVec(k);
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_A(A_rm, k); it_A; ++it_A)
            AB_ptr->insertBack(k, it_A.col()) = it_A.value();
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_B(B_rm, k); it_B; ++it_B)
            AB_ptr->insertBack(k, it_B.col() + A_rm.cols()) = it_B.value();
    }
    AB_ptr->finalize();
}