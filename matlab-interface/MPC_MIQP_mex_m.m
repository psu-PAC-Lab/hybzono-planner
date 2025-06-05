classdef MPC_MIQP_mex_m < cppclass

    properties (Access = private)

        % mi settings fields
        mi_settings_field_names = {...
            'eps_feas', ...
            'max_iter_bb', ...
            'verbose', ...
            'conv_rel', ...
            'conv_abs', ...
            'T_max', ...
            'n_threads'};

        % qp settings fields
        qp_settings_field_names = {...
            'mu_term', ...
            'mu_feas', ...
            'mu_max', ...
            'mu_init', ...
            'iter_max', ...
            'gamma', ...
            't_ls', ...
            'T_max'};

        % default mi settings
        mi_settings = struct(...
            'eps_feas', 1e-3, ...
            'max_iter_bb', 1e5, ...
            'verbose', true, ...
            'conv_rel', 1e-2, ...
            'conv_abs', 1e-1, ...
            'T_max', 0, ...
            'n_threads', 1);

        % default qp settings
        qp_settings = struct(...
            'mu_term', 1e-3, ... % convergence threshold
            'mu_feas', 1e0, ... % flag as feasible if mu < mu_feas
            'mu_max', 1e9, ... % max allowed value of barrier parameter during search (infeasibility detection)
            'mu_init', 1e5, ... % initial value of barrier parameter
            'iter_max', 100, ...
            'gamma', 0.999, ... % step size adjustment, <1
            't_ls', 0.9, ... % line search parameter in (0,1)
            'T_max', 0); 

        % mpc settings fields  
        mpc_settings_field_names = {...
            'warm_start', ...
            'n_horizon', ...
            'u1_control'};

        % default mpc settings
        mpc_settings = struct(...
            'warm_start', false, ...
            'n_horizon', 10, ...
            'u1_control', false);

    end

    % class interface
    methods (Access = public)

        function obj = MPC_MIQP_mex_m(varargin)
            obj@cppclass('MPC_MIQP_mex', varargin{:});
        end

        function set_dyn_matrices(obj, A_dyn, B_dyn)
            obj.cppmethod('set_dyn_matrices', sparse(A_dyn), sparse(B_dyn));
        end

        function set_dyn_matrices_LTV(obj, A_dyn_vec, B_dyn_vec)
            A_dyn_cell = cell(size(A_dyn_vec,3), 1);
            B_dyn_cell = cell(size(B_dyn_vec,3), 1);
            for ii = 1:size(A_dyn_vec,3)
                A_dyn_cell{ii} = sparse(A_dyn_vec(:,:,ii));
                B_dyn_cell{ii} = sparse(B_dyn_vec(:,:,ii));
            end

            obj.cppmethod('set_dyn_matrices_LTV', A_dyn_cell, B_dyn_cell);
        end

        function set_stage_cost(obj, Q, R, N, q_x, q_u)
            Q = sparse(Q);
            R = sparse(R);
            
            if nargin < 4
                N = [];
            else
                N = sparse(N);
            end
            if nargin < 5
                q_x = [];
            end
            if nargin < 6
                q_u = [];
            end

            obj.cppmethod('set_stage_cost', Q, R, N, q_x, q_u);
        end 

        function set_terminal_cost(obj, Qf, q_xN)
            Qf = sparse(Qf);
            if nargin < 3
                q_xN = [];
            end

            obj.cppmethod('set_terminal_cost', Qf, q_xN);
        end

        function set_state_inequality_constraints(obj, A_x, b_x)
            obj.cppmethod('set_state_inequality_constraints', sparse(A_x), b_x);
        end

        function set_input_inequality_constraints(obj, A_u, b_u)
            obj.cppmethod('set_input_inequality_constraints', sparse(A_u), b_u);
        end

        function set_terminal_state_inequality_constraints(obj, A_xf, b_xf)
            obj.cppmethod('set_terminal_state_inequality_constraints', sparse(A_xf), b_xf);
        end

        function set_state_inequality_constraints_conzono(obj, Zc_x, Pp_x, x_min_hard, x_max_hard)
            obj.cppmethod('set_state_inequality_constraints_conzono', obj.conzono2struct(Zc_x), full(Pp_x), x_min_hard, x_max_hard);
        end

        function set_input_inequality_constraints_conzono(obj, Zc_u, Pp_u, u_min_hard, u_max_hard)
            obj.cppmethod('set_input_inequality_constraints_conzono', obj.conzono2struct(Zc_u), full(Pp_u), u_min_hard, u_max_hard);
        end

        function set_terminal_state_inequality_constraints_conzono(obj, Zc_xf, Pp_xf, xf_min_hard, xf_max_hard)
            obj.cppmethod('set_terminal_state_inequality_constraints_conzono', obj.conzono2struct(Zc_xf), full(Pp_xf), xf_min_hard, xf_max_hard);
        end

        function set_state_inequality_constraints_slack_vars_cost(obj, Q_x, sigma_max_x)
            obj.cppmethod('set_state_inequality_constraints_slack_vars_cost', sparse(Q_x), sigma_max_x);
        end

        function set_input_inequality_constraints_slack_vars_cost(obj, Q_u, sigma_max_u)
            obj.cppmethod('set_input_inequality_constraints_slack_vars_cost', sparse(Q_u), sigma_max_u);
        end

        function set_terminal_state_inequality_constraints_slack_vars_cost(obj, Q_xf, sigma_max_xf)
            obj.cppmethod('set_terminal_state_inequality_constraints_slack_vars_cost', sparse(Q_xf), sigma_max_xf);
        end

        function set_mpc_horizon(obj, N)
            obj.cppmethod('set_mpc_horizon', N);
        end

        function set_qp_settings(obj, settings_in)
            obj.updateValidSettings_QP(settings_in);
            obj.cppmethod('set_qp_settings', obj.qp_settings);
        end

        function set_mi_settings(obj, settings_in)
            obj.updateValidSettings_MI(settings_in);
            obj.cppmethod('set_mi_settings', obj.mi_settings);
        end

        function set_mpc_settings(obj, settings_in)
            obj.updateValidSettings_MPC(settings_in);
            obj.cppmethod('set_mpc_settings', obj.mpc_settings);
        end

        function set_distance_limit(obj, d_max)
            obj.cppmethod('set_distance_limit', d_max);
        end

        function setup_hybzono_constraints(obj, varargin)
            if length(varargin) == 2
                X_hz = obj.hybzono2struct(varargin{1});
                Pp = varargin{2};

                obj.cppmethod('setup_hybzono_constraints', X_hz, Pp);

            elseif length(varargin) == 3
                X_hz = obj.hybzono2struct(varargin{1});
                Pp = varargin{2};
                Q_hz = sparse(varargin{3});

                obj.cppmethod('setup_hybzono_constraints', X_hz, Pp, Q_hz);

            elseif length(varargin) == 5
                X_hz = obj.hybzono2struct(varargin{1});
                Pp = varargin{2};
                Q_hz = sparse(varargin{3});
                sigma_hz_max = varargin{4};
                sigma_hz_min = varargin{5};

                obj.cppmethod('setup_hybzono_constraints', X_hz, Pp, Q_hz, sigma_hz_max, sigma_hz_min);

            else
                error('invalid number of input arguments');
            end
        end

        function setup_hybzono_Hrep_constraints(obj, X_hz, A_Hrep, b_Hrep, Pc, Q_hz, sigma_hz_min, sigma_hz_max, cost_vec, ...
                Pp, A_Hrep_pos, b_Hrep_pos)

            % mandatory arguments
            if isempty(X_hz) || isempty(A_Hrep) || isempty(b_Hrep) || isempty(Pc)
                error('missing mandatory arguments');
            end

            % optional arguments
            if nargin < 6
                Q_hz = [];
            end
            if nargin < 7
                sigma_hz_min = [];
            end
            if nargin < 8
                sigma_hz_max = [];
            end
            if nargin < 9
                cost_vec = [];
            end
            if nargin < 10
                Pp = [];
            end
            if nargin < 11
                A_Hrep_pos = [];
            end
            if nargin < 12
                b_Hrep_pos = [];
            end

            % process hybzono
            X_hz_st = obj.hybzono2struct(X_hz);

            % process H-rep
            A_Hrep_st.A_Hrep_mat = A_Hrep.get_mat_vec();
            A_Hrep_st.A_Hrep_ind = A_Hrep.get_ind_vec();

            b_Hrep_st.b_Hrep_mat = b_Hrep.get_mat_vec();
            b_Hrep_st.b_Hrep_ind = b_Hrep.get_ind_vec();

            if (~isempty(A_Hrep_pos))
                A_Hrep_pos_st.A_Hrep_pos_mat = A_Hrep_pos.get_mat_vec();
                A_Hrep_pos_st.A_Hrep_pos_ind = A_Hrep_pos.get_ind_vec();
            else
                A_Hrep_pos_st = [];
            end
            if (~isempty(b_Hrep_pos))
                b_Hrep_pos_st.b_Hrep_pos_mat = b_Hrep_pos.get_mat_vec();
                b_Hrep_pos_st.b_Hrep_pos_ind = b_Hrep_pos.get_ind_vec();
            else
                b_Hrep_pos_st = [];
            end
            
            % call cpp method
            obj.cppmethod('setup_hybzono_Hrep_constraints', X_hz_st, A_Hrep_st, b_Hrep_st, Pc, sparse(Q_hz), ...
                sigma_hz_min, sigma_hz_max, cost_vec, Pp, A_Hrep_pos_st, b_Hrep_pos_st);
            
        end

        function set_region_dependent_disturbances(obj, d_region)
            d_region_struct.d_region_vec = d_region.get_mat_vec();
            d_region_struct.d_region_ind = d_region.get_ind_vec();
            obj.cppmethod('set_region_dependent_disturbances', d_region_struct);
        end

        function set_dist_between_regions(obj, dist_mat)
            obj.cppmethod('set_dist_between_regions', dist_mat);
        end

        function set_terminal_hybzono_constraints(obj, X_hz_term, term_region_cost_vec, Pp_term, Q_term)
            obj.cppmethod('set_terminal_hybzono_constraints', obj.hybzono2struct(X_hz_term), term_region_cost_vec, Pp_term, sparse(Q_term));
        end

        function build_controller(obj)
            obj.cppmethod('build_controller');
        end

        function varargout = control(obj, x, x_ref)
            [u, feasible] = obj.cppmethod('control', x, x_ref);
            outputs = {u, feasible};
            for ii = 1:nargout
                varargout{ii} = outputs{ii};
            end
        end

        function varargout = control_LTV(obj, x, x_ref, A_dyn_vec, B_dyn_vec)
            A_dyn_cell = cell(size(A_dyn_vec,3), 1);
            B_dyn_cell = cell(size(B_dyn_vec,3), 1);
            for ii = 1:size(A_dyn_vec,3)
                A_dyn_cell{ii} = sparse(A_dyn_vec(:,:,ii));
                B_dyn_cell{ii} = sparse(B_dyn_vec(:,:,ii));
            end
            [u, feasible] = obj.cppmethod('control_LTV', x, x_ref, A_dyn_cell, B_dyn_cell);
            outputs = {u, feasible};
            for ii = 1:nargout
                varargout{ii} = outputs{ii};
            end
        end

        function [x_vec, u_vec] = get_trajectory(obj)
            [x_vec, u_vec] = obj.cppmethod('get_trajectory');
        end
        
        function region = get_terminal_region_selection(obj)
            region = obj.cppmethod('get_terminal_region_selection');
        end

        function varargout = get_miqp_results(obj)
            [results, x, region_vec] = obj.cppmethod('get_miqp_results');
            outputs = {results, x, region_vec};
            for ii = 1:nargout
                varargout{ii} = outputs{ii};
            end
        end

        function do_u1_control(obj, u0, x0)
            obj.cppmethod('do_u1_control', u0, x0);
        end

    end



    % private methods
    methods (Access = private)

        % hybzono to struct
        function X_hz_struct = hybzono2struct(obj, X_hz)

            X_hz_struct.Gc = X_hz.Gc;
            X_hz_struct.Gb = X_hz.Gb;
            X_hz_struct.c = X_hz.c;
            X_hz_struct.Ac = X_hz.Ac;
            X_hz_struct.Ab = X_hz.Ab;
            X_hz_struct.b = X_hz.b;
            X_hz_struct.n = X_hz.n;
            X_hz_struct.nGc = X_hz.nGc;
            X_hz_struct.nGb = X_hz.nGb;

        end

        % conzono to struct
        function Zc_struct = conzono2struct(obj, Zc)
            
            Zc_struct.G = Zc.G;
            Zc_struct.c = Zc.c;
            Zc_struct.A = Zc.A;
            Zc_struct.b = Zc.b;
            Zc_struct.n = Zc.n;
            Zc_struct.nG = Zc.nG;
            Zc_struct.nC = Zc.nC;

        end

        % update valid MI settings
        function updateValidSettings_MI(obj, settings_in)
            
            % get fieldnames
            fldnames = fieldnames(settings_in);

            % loop through fields
            for ii = 1:length(fldnames)
                
                % check if field name is valid
                if any(strcmp(obj.mi_settings_field_names, fldnames{ii}))
                    
                    % update field in settings structure
                    new_val = getfield(settings_in, fldnames{ii});
                    obj.mi_settings = setfield(obj.mi_settings, fldnames{ii}, new_val);

                end
            end
        end

        % update valid QP settings
        function updateValidSettings_QP(obj, settings_in)
            
            % get fieldnames
            fldnames = fieldnames(settings_in);

            % loop through fields
            for ii = 1:length(fldnames)
                
                % check if field name is valid
                if any(strcmp(obj.qp_settings_field_names, fldnames{ii}))
                    
                    % update field in settings structure
                    new_val = getfield(settings_in, fldnames{ii});
                    obj.qp_settings = setfield(obj.qp_settings, fldnames{ii}, new_val);

                end
            end
        end

        % update valid MPC settings
        function updateValidSettings_MPC(obj, settings_in)

            % get fieldnames
            fldnames = fieldnames(settings_in);

            % loop through fields
            for ii = 1:length(fldnames)

                % check if field name is valid
                if any(strcmp(obj.mpc_settings_field_names, fldnames{ii}))

                    % update field in settings structure
                    new_val = getfield(settings_in, fldnames{ii});
                    obj.mpc_settings = setfield(obj.mpc_settings, fldnames{ii}, new_val);

                end
            end
        end


    end


end