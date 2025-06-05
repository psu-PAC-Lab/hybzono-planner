#include "MI_MPC.hpp"

using namespace MI_QPIPMPC;

// move constructor
MI_MPC::MI_MPC(MI_MPC &&other)
{
    update_pending = other.update_pending;
    data = std::move(other.data);
    data_volatile = std::move(other.data_volatile);
    verbose_output = std::move(other.verbose_output);
}

// copy assignment
MI_MPC& MI_MPC::operator=(const MI_MPC &other)
{
    if (this != &other)
    {
        update_pending = other.update_pending;
        data = std::move(other.data);
        data_volatile = std::move(other.data_volatile);
        verbose_output = std::move(other.verbose_output);
    }
    return *this;
}

// setup method
void MI_MPC::setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
    const QP_IP_MPC::QP_settings &qp_settings, 
    const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
    const std::map<int, std::vector<int>> &R_x02region, int region_0,
    const std::map<int, Eigen::MatrixXd> &A_Hrep,
    const std::map<int, Eigen::VectorXd> &b_Hrep)
{
    // start setup timer
    auto start = std::chrono::high_resolution_clock::now();

    // setup data
    data = Data(qpipmpc_data, R_x02region, R_region2region, region_0, settings, qp_settings, A_Hrep, b_Hrep);
    data_volatile = DataVolatile(data);

    // record setup time
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_final = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    data_volatile.set_setup_time(1e-6 * ((double) duration_final.count()));
}

// setup method
void MI_MPC::setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
    const QP_IP_MPC::QP_settings &qp_settings, 
    const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
    const std::map<int, std::vector<int>> &R_x02region, int region_0)
{
    // start setup timer
    auto start = std::chrono::high_resolution_clock::now();

    // setup data
    data = Data(qpipmpc_data, R_x02region, R_region2region, region_0, settings, qp_settings);
    data_volatile = DataVolatile(data);

    // record setup time
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_final = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    data_volatile.set_setup_time(1e-6 * ((double) duration_final.count()));
}

void MI_MPC::setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            const std::map<int, std::vector<int>> &R_x02region, int region_0,
            const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep,
            const Eigen::SparseMatrix<double> &C_term,
            const Eigen::SparseMatrix<double> &D_term,
            const Eigen::SparseMatrix<double> &G_term,
            const Eigen::SparseMatrix<double> &P_term,
            const Eigen::Ref<const Eigen::VectorXd> crhs_term,
            const Eigen::Ref<const Eigen::VectorXd> w_term,
            const Eigen::Ref<const Eigen::VectorXd> q_term,
            const std::vector<int> &idx_term,
            const std::vector<int> &idx_binvar_term)
{
    // start setup timer
    auto start = std::chrono::high_resolution_clock::now();

    // setup data
    data = Data(qpipmpc_data, R_x02region, R_region2region, region_0, settings, qp_settings, A_Hrep, b_Hrep);
    data_volatile = DataVolatile(data);

    // set terminal constraint
    data.set_terminal_constraint(C_term, D_term, G_term, P_term, crhs_term, w_term, q_term, idx_term, idx_binvar_term);

    // record setup time
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_final = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    data_volatile.set_setup_time(1e-6 * ((double) duration_final.count()));
}

void MI_MPC::setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            const std::map<int, std::vector<int>> &R_x02region, int region_0,
            const Eigen::SparseMatrix<double> &C_term,
            const Eigen::SparseMatrix<double> &D_term,
            const Eigen::SparseMatrix<double> &G_term,
            const Eigen::SparseMatrix<double> &P_term,
            const Eigen::Ref<const Eigen::VectorXd> crhs_term,
            const Eigen::Ref<const Eigen::VectorXd> w_term,
            const Eigen::Ref<const Eigen::VectorXd> q_term,
            const std::vector<int> &idx_term,
            const std::vector<int> &idx_binvar_term)
{
    // start setup timer
    auto start = std::chrono::high_resolution_clock::now();

    // setup data
    data = Data(qpipmpc_data, R_x02region, R_region2region, region_0, settings, qp_settings);
    data_volatile = DataVolatile(data);

    // set terminal constraint
    data.set_terminal_constraint(C_term, D_term, G_term, P_term, crhs_term, w_term, q_term, idx_term, idx_binvar_term);

    // record setup time
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_final = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    data_volatile.set_setup_time(1e-6 * ((double) duration_final.count()));
}

// solve method
Results MI_MPC::solve()
{
    using namespace std::chrono_literals;

    // start solve timer
    auto start = std::chrono::high_resolution_clock::now();

    // update if necessary
    if (update_pending)
    {
        update();
        update_pending = false;
    }

    // verbosity
    if (data.settings.verbose)
        verbose_output.print_headline();
    
    // time limit for optimization
    double T_max;
    bool enforce_max_time;
    if (data.settings.T_max > 0)
    {
        T_max = data.settings.T_max;
        enforce_max_time = true;
    }
    else
    {
        T_max = QP::inf;
        enforce_max_time = false;
    }
    double running_timer = 0;

    // get maximum number of threads
    // -1 on hardware concurrency to account for the thread running this process
    const int n_threads = std::min(data.settings.n_threads, std::thread::hardware_concurrency()-1);

    // init vector holding the lower bounds of the nodes in the active threads
    data_volatile.init_lower_glob_threads(n_threads);

    if (n_threads > 1)
    {
        // init vectors of threads, futures, real_threads flag
        std::vector<std::tuple<std::thread, std::future<void>, bool>> thread_data;
        thread_data.reserve(n_threads);

        // start threads
        std::packaged_task<void(int)> task; // declare
        std::future<void> future; // declare
        int n_real_threads = std::min(n_threads, leaves.size());
        for (int i=0; i<n_real_threads; i++)
        {
            task = std::packaged_task<void(int)>([this] (int i) {main_solve_loop(i);});
            future = task.get_future();
            thread_data.emplace_back(std::make_tuple(std::thread(std::move(task), i), std::move(future), true));
        }
        for (int i=n_real_threads; i<n_threads; i++)
        {
            thread_data.emplace_back(std::make_tuple(std::thread(), std::future<void>(), false));
        }

        int n_dead_threads = n_threads - n_real_threads;
        bool can_start_new_thread = false; // init

        // solve: loop until there are no more leaves, converged, or hit max time
        do // must run at least once
        {   
            // init number of dead threads
            n_dead_threads = 0;

            // loop over all threads
            auto it = thread_data.begin();
            for (int i=0; i<n_threads; i++)
            {
                // check if thread complete
                if (std::get<0>(*it).joinable()) // if thread is joinable
                {
                    auto status = std::get<1>(*it).wait_for(0ms); // get future status
                    if (status == std::future_status::ready)
                    {
                        // join thread
                        std::get<0>(*it).join();

                        // reset lower_glob_threads
                        if (std::get<2>(*it))
                            data_volatile.reset_lower_glob_threads(i);

                        // mark thread as able to be restarted
                        can_start_new_thread = true;
                    }
                    else
                    {
                        can_start_new_thread = false;
                    }
                }
                else
                {
                    can_start_new_thread = true;
                }

                if (can_start_new_thread)
                {
                    // if there are leaves, create new task, future, and thread
                    if (!leaves.empty())
                    {
                        task = std::packaged_task<void(int)>([this] (int i) {main_solve_loop(i);});
                        std::get<1>(*it) = task.get_future();
                        std::get<0>(*it) = std::thread(std::move(task), i);
                        std::get<2>(*it) = true;
                    }    
                    else
                    {
                        std::get<2>(*it) = false;
                    }
                }

                // count number of dead threads
                if (!(std::get<2>(*it)))
                    ++n_dead_threads;

                // increment
                ++it;
            }    

            // running timer
            if (enforce_max_time)
            {
                auto timer_running = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timer_running - start);
                running_timer = 1e-6 * ((double) duration.count());
            }
        }
        while (can_continue(n_dead_threads, n_threads) && !is_converged() && (running_timer < T_max));

        // join any remaining threads
        for (auto it = thread_data.begin(); it != thread_data.end(); ++it)
        {
            if (std::get<0>(*it).joinable())
                std::get<0>(*it).join();
        }
    }
    else // single-threaded
    {
        // solve: loop until there are no more leaves, converged, or hit max time
        do
        {
            main_solve_loop(0);

            // running timer
            if (enforce_max_time)
            {
                auto timer_running = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timer_running - start);
                running_timer = 1e-6 * ((double) duration.count());
            }
        }
        while (data_volatile.get_iter_num() < data.settings.max_iter_bb && !leaves.empty() && !is_converged() && (running_timer < T_max));
    }

    // get final status
    update_return_status();

    // update average number of iterations
    data_volatile.calc_qp_iter_avg();

    // get solve time
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_final = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    data_volatile.set_solve_time(1e-6 * ((double) duration_final.count()));

    // calculate run time
    data_volatile.calc_run_time();

    // print footer
    if (data.settings.verbose)
        verbose_output.print_footer(data_volatile);

    // return results
    return Results(data_volatile.get_x(), data_volatile.get_upper_glob(), data_volatile.get_lower_glob(),
        data_volatile.get_run_time(), data_volatile.get_status(), data_volatile.get_qp_solve_time(),
        data_volatile.get_qp_iter_avg(), data_volatile.get_iter_num(), data_volatile.get_region_vec());
}

// main solver loop
void MI_MPC::main_solve_loop(int thread_number)
{
    // 1) choose leaf
    std::unique_ptr<Node> leaf = leaves.pop_top(data_volatile, thread_number);

    // check if valid
    if (!leaf->is_valid()) return;

    // 2) solve relaxed problem in leaf
    leaf->solve(data);

    // 3) branch and bound
    branch_and_bound(leaf);

    // verbosity
    if (data.settings.verbose)
        verbose_output.print_progress(leaf, leaves, data_volatile);

    // update iteration number
    data_volatile.increment_iter_num();
}


// core branch-and-bound function
void MI_MPC::branch_and_bound(const std::unique_ptr<Node> &leaf)
{
    // update total number of iterations
    data_volatile.increment_qp_iter(leaf->num_iter);

    // update total time to solve QP problems
    data_volatile.increment_qp_solve_time(leaf->qp_solve_time);

    // 1) if infeasible or unbounded, then return (prune)
    if (leaf->status == QP_INFEASIBLE)
        return;
    
    // 2) if lower bound greater than upper bound, then return (prune)
    if (leaf->lower > data_volatile.get_upper_glob())
        return;

    // 3) if leaf is full-depth, store its solution if optimal
    bool full_depth; // declare
    if (data.terminal_constraint)
        full_depth = (leaf->depth == data.n_horizon+1) && (leaf->term_constraint_depth == 1);
    else
        full_depth = (leaf->depth == data.n_horizon+1);

    if (full_depth)
    {
        // check whether leaf beats current optimum
        if (leaf->obj < data_volatile.get_upper_glob())
        {
            // update optimum
            data_volatile.update_optimum(leaf);

            // prune
            leaves.prune(data_volatile.get_upper_glob());
        }

        return;
    }

    // 4) if not full depth, no collisions detected, and estimated objective less than incumbent, generate leaf and add to top of queue
    bool can_do_collision_checking, leaf_valid; // declare
    double approx_cost; // declare 
    std::unique_ptr<Node> copy_leaf (new Node(*leaf)); // init

    if (data.Hrep_defined)
    {
        if (data.terminal_constraint)
            can_do_collision_checking = (leaf->depth < data.n_horizon+1) && (leaf->term_constraint_depth == 1);
        else
            can_do_collision_checking = (leaf->depth < data.n_horizon+1);
    }
    else
    {
        can_do_collision_checking = false;
    }

    if (can_do_collision_checking)
    {
        leaf->get_regions_no_constraint_violation(data);
        if (!leaf->cons_violated)
        {
            // get approximate cost assicated with choosing regions in region_vec_no_violattion
            approx_cost = leaf->get_approx_cost(data, leaf->region_vec_no_violation);

            // if cost less than incumbent, create leaf and add to front of queue
            if (approx_cost < data_volatile.get_upper_glob())
            {
                // copy region_vec_no_violation into region_vec
                copy_leaf->region_vec = leaf->region_vec_no_violation;

                // update leaf
                copy_leaf->apply_reachability_constraints(data);
                leaf_valid = copy_leaf->check_reachability_valid(data);
                if (leaf_valid)
                {
                    copy_leaf->get_depth(data);
                    copy_leaf->give_priority();
                    
                    // add leaf to front of queue
                    leaves.push_top(std::move(copy_leaf));
                }

                // do not return, continue to branch
            }
        }
    }

    // 5) if we got here, branch
    if (data.terminal_constraint && leaf->term_constraint_depth == 0)
    {   
        // branch at terminal constraint first
        branch_at_terminal_constraint(leaf);
    }
    else if (!data.Hrep_defined)
    {
        throw std::runtime_error("branching without H-rep not yet supported");
    }
    else if (leaf->cons_violated)
    {
        branch_at_timestep(leaf, leaf->k_first_cons_violation);
    }
    else
    {
        branch_most_fractional_Hrep(leaf);
    }
        
    // 6) update global lower bound
    data_volatile.set_lower_glob(leaves.get_lower_glob());
}


// branching logic: branch at obstacle violation
void MI_MPC::branch_at_timestep(const std::unique_ptr<Node> &leaf, int kb) 
{
    // if leaf is full depth, return
    if (leaf->depth == data.n_horizon+1)
        return;
    
    // get fractionality of each region in solution at kb
    Eigen::VectorXd frac_vec = leaf->get_frac_vec(kb, leaf->region_vec.at(kb), data);

    // get branching region
    int region_b = -1; // init
    double frac_max = frac_vec.maxCoeff();
    int cnt = 0;
    while ((region_b == -1) && (cnt < frac_vec.size()))
    {
        if (frac_vec[cnt] == frac_max)
            region_b = leaf->region_vec.at(kb)[cnt];
        cnt++;
    }

    // create branch leaf leaf_b
    std::unique_ptr<Node> leaf_b (new Node(*leaf));
    bool leaf_b_valid = leaf_b->update_reachability_constraints(data, kb, {region_b});
    if (leaf_b_valid)
    {
        leaf_b->get_depth(data);
        leaf_b->set_frac(frac_max);
        leaf_b->remove_priority();
    }

    // get remaining regions
    std::vector<int> regions_r = leaf->region_vec.at(kb);
    regions_r.erase(std::remove(regions_r.begin(), regions_r.end(), region_b), 
        regions_r.end());

    // create remeaining leaf leaf_r
    std::unique_ptr<Node> leaf_r (new Node(*leaf));
    bool leaf_r_valid = leaf_r->update_reachability_constraints(data, kb, regions_r);
    if (leaf_r_valid)
    {
        leaf_r->get_depth(data);
        Eigen::VectorXd frac_vec_r = leaf_r->get_frac_vec(kb, regions_r, data);
        leaf_r->set_frac(frac_vec_r.maxCoeff());
        leaf_r->remove_priority();
    }

    // check that upper bound has not been updated and add leaves if applicable
    if (data_volatile.get_upper_glob() > leaf->lower)
    {
        // add leaves to queue
        if (leaf_b_valid)
            leaves.push(std::move(leaf_b));
        if (leaf_r_valid)
            leaves.push(std::move(leaf_r));
    }
}

// branching logic: branch at terminal constraint
void MI_MPC::branch_at_terminal_constraint(const std::unique_ptr<Node> &leaf)
{
    // if leaf is full depth (wrt terminal constraint), return
    if (leaf->term_constraint_depth == 1)
        return;

    // get fractionality of each terminal constraint region in solution
    Eigen::VectorXd frac_vec = leaf->get_term_frac_vec(leaf->term_region_vec, data);

    // get branching region
    int region_b = -1; // init
    double frac_max = frac_vec.maxCoeff();
    int cnt = 0;
    while ((region_b == -1) && (cnt < frac_vec.size()))
    {
        if (frac_vec[cnt] == frac_max)
            region_b = leaf->term_region_vec[cnt];
        cnt++;
    }

    // create branch leaf leaf_b
    std::unique_ptr<Node> leaf_b (new Node(*leaf));
    leaf_b->update_terminal_constraint({region_b});
    leaf_b->get_terminal_constraint_depth();
    leaf_b->set_term_frac(frac_max);
    leaf_b->remove_priority();

    // get remaining regions
    std::vector<int> regions_r = leaf->term_region_vec;
    regions_r.erase(std::remove(regions_r.begin(), regions_r.end(), region_b), 
        regions_r.end());

    // create remeaining leaf leaf_r
    std::unique_ptr<Node> leaf_r (new Node(*leaf));
    leaf_r->update_terminal_constraint(regions_r);
    leaf_r->get_terminal_constraint_depth();
    Eigen::VectorXd frac_vec_r = leaf_r->get_term_frac_vec(regions_r, data);
    leaf_r->set_term_frac(frac_vec_r.maxCoeff());
    leaf_r->remove_priority();

    // check that upper bound has not been updated and add leaves if applicable
    if (data_volatile.get_upper_glob() > leaf->lower)
    {
        // add leaves to queue
        leaves.push(std::move(leaf_b));
        leaves.push(std::move(leaf_r));
    }
}

void MI_MPC::branch_most_fractional_Hrep(const std::unique_ptr<Node> &leaf)
{
    // if leaf is full depth, return
    if (leaf->depth == data.n_horizon+1)
        return;
    
    // do most fractional branching

    // declare/init working vars
    Eigen::VectorXd frac_vec;
    double max_dist_from_int = 0;
    int region_b = -1;
    std::vector<int> regions_r;
    regions_r.reserve(data.n_regions);
    double dist_from_int, frac_max;
    int kb = -1;
    int cnt;

    // find time step and region with most fractional value
    for (int k=0; k<=data.n_horizon; ++k)
    {
        // fractionality of each region at time step k
        frac_vec = leaf->get_frac_vec(k, leaf->region_vec.at(k), data);
        frac_max = frac_vec.maxCoeff();
        dist_from_int = std::abs(frac_max - std::round(frac_max));

        if (dist_from_int > max_dist_from_int)
        {
            // update incumbent branching timestep and regions
            max_dist_from_int = dist_from_int;
            kb = k;
        }
    }
 
    // branch at timestep
    branch_at_timestep(leaf, kb);
}

// update methods
void MI_MPC::update()
{
    // clear leaves queue
    leaves.clear();

    // declare terminal region vector
    std::vector<int> term_region_vec; 
    term_region_vec.reserve(data.n_term_regions);

    // warm-start
    int n_warm_start_nodes;
    std::map<int, std::vector<int>> region_vec_warm_start;
    std::unique_ptr<Node> warm_start_node; // declare
    if (data.warm_start)
    {
        // create warm start nodes
        n_warm_start_nodes = data.region_vec_warm_start.size();

        for (int i=0; i<n_warm_start_nodes; i++)
        {
            // get region_vec
            region_vec_warm_start.clear();
            for (int k=0; k<=data.n_horizon; k++)
            {
                region_vec_warm_start[k].push_back(data.region_vec_warm_start[i][k]);
            }

            // make node
            if (data.terminal_constraint)
            {
                term_region_vec.clear();
                term_region_vec.push_back(data.term_region_warm_start);

                warm_start_node = std::make_unique<Node>(data, data.solver, 
                    -QP::inf, region_vec_warm_start, term_region_vec);
            }
            else
            {
                warm_start_node = std::make_unique<Node>(data, data.solver, 
                    -QP::inf, region_vec_warm_start);

            }
            warm_start_node->reset_solver(data);
            warm_start_node->apply_reachability_constraints(data);
            warm_start_node->get_depth(data);
            warm_start_node->give_priority();

            // add to leaves queue
            if (warm_start_node->check_reachability_valid(data))
            {
                leaves.push(std::move(warm_start_node));
            }
        }
    }

    // create root node
    std::unique_ptr<Node> root(new Node); // allocate
    std::map<int, std::vector<int>> region_vec_empty;
    if (data.terminal_constraint)
    {
        term_region_vec.clear();
        term_region_vec.resize(data.n_term_regions);
        std::iota(term_region_vec.begin(), term_region_vec.end(), 1); // set to 1, 2, ..., n_term_regions

        *root = Node(data, data.solver, 
            -QP::inf, region_vec_empty, term_region_vec);
    }
    else
    {
        *root = Node(data, data.solver, 
            -QP::inf, region_vec_empty);
    }
    root->reset_solver(data);
    bool root_valid = root->init_reachability_constraints(data);
    root->give_priority();
    
    // add root node to leaves - make sure ahead of any warm start nodes
    if (root_valid)
        leaves.push_top(std::move(root));

    // reset problem data
    data_volatile.reset(data);

    // reset verbose output
    if (data.settings.verbose)
        verbose_output.reset();
}

void MI_MPC::update_P_i_nom(const Eigen::SparseMatrix<double> &P_i_nom)
{
    data.set_P_i_nom(P_i_nom);
    update_pending = true;
}

void MI_MPC::update_q_i_nom(const Eigen::Ref<const Eigen::VectorXd> q_i_nom)
{
    data.set_q_i_nom(q_i_nom);
    update_pending = true;
}

void MI_MPC::update_b(double b)
{
    data.set_b(b);
    update_pending = true;
}

void MI_MPC::update_C_i_nom(const Eigen::SparseMatrix<double> &C_i_nom)
{
    data.set_C_i_nom(C_i_nom);
    update_pending = true;
}

void MI_MPC::update_D_i_nom(const Eigen::SparseMatrix<double> &D_i_nom)
{
    data.set_D_i_nom(D_i_nom);
    update_pending = true;
}

void MI_MPC::update_crhs_i_nom(const Eigen::Ref<const Eigen::VectorXd> crhs_i_nom)
{
    data.set_crhs_i_nom(crhs_i_nom);
    update_pending = true;
}

void MI_MPC::update_G_i_nom(const Eigen::SparseMatrix<double> &G_i_nom)
{
    data.set_G_i_nom(G_i_nom);
    update_pending = true;
}

void MI_MPC::update_w_i_nom(const Eigen::Ref<const Eigen::VectorXd> w_i_nom)
{
    data.set_w_i_nom(w_i_nom);
    update_pending = true;
}

void MI_MPC::update_P_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &P_i_vec)
{
    data.set_P_i_vec(P_i_vec);
    update_pending = true;
}

void MI_MPC::update_q_i_vec(const std::map<int, Eigen::VectorXd> &q_i_vec)
{
    data.set_q_i_vec(q_i_vec);
    update_pending = true;
}

void MI_MPC::update_C_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &C_i_vec)
{
    data.set_C_i_vec(C_i_vec);
    update_pending = true;
}

void MI_MPC::update_D_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &D_i_vec)
{
    data.set_D_i_vec(D_i_vec);
    update_pending = true;
}

void MI_MPC::update_crhs_i_vec(const std::map<int, Eigen::VectorXd> &crhs_i_vec)
{
    data.set_crhs_i_vec(crhs_i_vec);
    update_pending = true;
}

void MI_MPC::update_G_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &G_i_vec)
{
    data.set_G_i_vec(G_i_vec);
    update_pending = true;
}

void MI_MPC::update_w_i_vec(const std::map<int, Eigen::VectorXd> &w_i_vec)
{
    data.set_w_i_vec(w_i_vec);
    update_pending = true;
}

void MI_MPC::update_idx_i_nom(const std::vector<int> &idx_i_nom)
{
    data.set_idx_i_nom(idx_i_nom);
    update_pending = true;
}

void MI_MPC::update_idx_i_vec(const std::map<int, std::vector<int>> &idx_i_vec)
{
    data.set_idx_i_vec(idx_i_vec);
    update_pending = true;
}

void MI_MPC::update_idx_state(const std::map<int, std::vector<int>> &idx_state)
{
    data.set_idx_state(idx_state);
    update_pending = true;
}

void MI_MPC::update_idx_input(const std::map<int, std::vector<int>> &idx_input)
{
    data.set_idx_input(idx_input);
    update_pending = true;
}

void MI_MPC::update_idx_binvar(const std::map<int, std::vector<int>> &idx_binvar)
{
    data.set_idx_binvar(idx_binvar);
    update_pending = true;
}

void MI_MPC::update_idx_y(const std::map<int, std::vector<int>> &idx_y)
{
    data.set_idx_y(idx_y);
    update_pending = true;
}

void MI_MPC::update_idx_eq(const std::map<int, std::vector<int>> &idx_eq)
{
    data.set_idx_eq(idx_eq);
    update_pending = true;
}

void MI_MPC::update_idx_ineq(const std::map<int, std::vector<int>> &idx_ineq)
{
    data.set_idx_ineq(idx_ineq);
    update_pending = true;
}

void MI_MPC::update_reachability(const std::map<int, std::vector<int>> &R_x02region,
    const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
    int region_0)
{
    data.update_reachability(R_x02region, R_region2region, region_0);
    update_pending = true;
}

void MI_MPC::update_Hrep(const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep)
{
    data.set_Hrep(A_Hrep, b_Hrep);
    update_pending = true;
}

void MI_MPC::update_terminal_constraint(const Eigen::SparseMatrix<double> &C_term,
    const Eigen::SparseMatrix<double> &D_term,
    const Eigen::SparseMatrix<double> &G_term,
    const Eigen::SparseMatrix<double> &P_term,
    const Eigen::Ref<const Eigen::VectorXd> crhs_term,
    const Eigen::Ref<const Eigen::VectorXd> w_term,
    const Eigen::Ref<const Eigen::VectorXd> q_term,
    const std::vector<int> &idx_term,
    const std::vector<int> &idx_binvar_term)
{
    data.set_terminal_constraint(C_term, D_term, G_term, P_term, crhs_term, w_term, q_term, idx_term, idx_binvar_term);
    update_pending = true;
}

// get logging string when verbose
std::string MI_MPC::get_log() const
{
    return verbose_output.get_output_string();
}

void MI_MPC::set_output_stream_fcn(void (*funcptr)(const std::string &str))
{
    verbose_output.set_output_stream_fcn(funcptr);
}

// update return status
void MI_MPC::update_return_status() 
{
    double upper_glob = data_volatile.get_upper_glob();
    double lower_glob = data_volatile.get_lower_glob();
    int iter_num = data_volatile.get_iter_num();

    if (lower_glob > data.settings.max_cost)
    {
        data_volatile.set_status(MIQP_MAX_COST);
    }
    else if (iter_num < data.settings.max_iter_bb && (leaves.empty() || is_converged())) // finished
    {
        if (upper_glob < QP::inf && upper_glob > -QP::inf)
            data_volatile.set_status(MIQP_SOLVED);
        else
            data_volatile.set_status(MIQP_INFEASIBLE);
    }
    else // hit max number of iterations
    {   if (upper_glob < QP::inf && upper_glob > -QP::inf)
            data_volatile.set_status(MIQP_MAX_ITER_FEASIBLE);
        else
            data_volatile.set_status(MIQP_MAX_ITER_INFEASIBLE);
    }
}

// check if converged
bool MI_MPC::is_converged() const 
{
    double upper_glob = data_volatile.get_upper_glob();
    double lower_glob = data_volatile.get_lower_glob();

    if ((std::abs(upper_glob - lower_glob)/std::abs(upper_glob) < data.settings.conv_rel) ||
        (std::abs(upper_glob - lower_glob) < data.settings.conv_abs) ||
        (lower_glob > data.settings.max_cost))
        return true;
    else
        return false;
}

// check if can continue
bool MI_MPC::can_continue(int n_dead_threads, int n_threads) const 
{
    int iter_num = data_volatile.get_iter_num();

    if (n_dead_threads == n_threads || iter_num >= data.settings.max_iter_bb)
        return false;
    else
        return true;
}

// set warm start
void MI_MPC::set_warm_start_regions(const std::vector<std::vector<int>> &region_vec_warm_start, int term_region_warm_start)
{
    data.set_warm_start_regions(region_vec_warm_start, term_region_warm_start);
}