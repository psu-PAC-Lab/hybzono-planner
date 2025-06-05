classdef Map < matlab.mixin.Copyable
    
    properties (Access = protected)
        
        ind_vec = [];
        mat_vec = {};

    end

    methods (Access = public)
        
        % constructor
        function obj = Map()
        end

        % assignment
        function obj = subsasgn(obj, ind, mat)
            % input processing: https://www.mathworks.com/help/fixedpoint/ref/embedded.fi.subsasgn.html
            ind = ind.subs{1};
            ind_srch = find(obj.ind_vec == ind, 1);
            if isempty(ind_srch)
                obj.mat_vec{length(obj.mat_vec)+1} = mat;
                if isempty(obj.ind_vec)
                    obj.ind_vec = ind;
                else
                    obj.ind_vec(end+1) = ind;
                end
            else
                obj.mat_vec{ind_srch} = mat;
            end
        end

        % get matrix and method handling
        function varargout = subsref(obj, ind)
            switch ind(1).type
                case '()'
                    ind_val = ind(1).subs{1};
                    mat = obj.getVal(ind_val);
                    varargout{1} = mat;
                case '.'
                    methodname = ind(1).subs;
                    switch methodname
                        case 'exists'   
                            val = ind(2).subs{1};
                            varargout{1} = obj.exists(val);
                        case 'copy'
                            varargout{1} = obj.copy();
                        case 'get_mat_vec'
                            varargout{1} = obj.get_mat_vec();
                        case 'get_ind_vec'
                            varargout{1} = obj.get_ind_vec();
                        case 'length'
                            varargout{1} = obj.length();
                        case 'clear'
                            val = ind(2).subs{1};
                            varargout{1} = obj.clear(val);
                        otherwise
                            error('unknown method')
                    end
                otherwise
                    error('subsref error')
            end
        end

        % check if value exists
        function val_exists = exists(obj, ind)
            % see if index is in ind_vec
            ind_srch = find(obj.ind_vec == ind, 1);
            val_exists = ~isempty(ind_srch);
        end

        % get value at index
        function val = getVal(obj, ind)
            ind_srch = find(obj.ind_vec == ind, 1);
            if ~isempty(ind_srch)
                val = obj.mat_vec{ind_srch};
            else
                val = [];
            end
        end

        % get data
        function mat_vec = get_mat_vec(obj)
            mat_vec = obj.mat_vec;
        end
        
        function ind_vec = get_ind_vec(obj)
            ind_vec = obj.ind_vec;
        end

        % overridden "length"
        function len = length(obj)
            len = length(obj.ind_vec);
        end

        % overridden equality check
        function out = eq(obj, other)
            if length(obj) == length(other)
                out = true; % init
                for ii = obj.get_ind_vec()
                    if ~other.exists(ii) || (length(obj.getVal(ii)) ~= length(other.getVal(ii))) || ...
                        any(any(obj.getVal(ii) ~= other.getVal(ii)))
                        out = false;
                        return;
                    end
                end
            else
                out = false;
            end
        end

        % clear entry
        function out = clear(obj, ind)
            ind_srch = find(obj.ind_vec == ind, 1);
            ind_new = setdiff(1:length(obj.ind_vec), ind_srch);
            obj.ind_vec = obj.ind_vec(ind_new);
            obj.mat_vec = {obj.mat_vec{ind_new}};
            out = [];
        end

    end


end