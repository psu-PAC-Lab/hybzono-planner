%% setup

clear all
close all

% numbers of timesteps
timesteps = {'n_horizon: 5', 'n_horizon: 10', 'n_horizon: 15', 'n_horizon: 20'};
timesteps_num = [5,10,15,20];
timesteps_num_str = {'5', '10', '15', '20'};

% map names
maps = {'poly2', 'nonconvpoly1', 'zonopsu', 'zonocost1', 'zono3D'};

% solver names
solvers = {'ours_hz_warmstart', 'ours_hz', 'ours_hrep', 'gurobi_hz', 'gurobi_hrep'};

% colors associated with solvers
face_colors{1} = [0,0.6,1];
face_colors{2} = [0,0,1];
face_colors{3} = [0,0,1];
face_colors{4} = [1,0,0];
face_colors{5} = [1,0,0];

face_alphas{1} = 1;
face_alphas{2} = 1;
face_alphas{3} = 0.4;
face_alphas{4} = 1;
face_alphas{5} = 0.4;

% map and solver names as they appears in the plots
maps_plot = {'\textbf{(a)}', '\textbf{(b)}', '\textbf{(c)}', '\textbf{(d)}', '\textbf{(e)}'};
solvers_plot = {'Ours HZ WS', 'Ours HZ', 'Ours H-rep', 'Gurobi HZ', 'Gurobi H-rep'};


%% read text file

fid_2D = fopen('solver_comparison.txt', 'r');
fid_3D = fopen('solver_comparison_3D.txt', 'r');
on_3D = false;

cnt_max = 1e4;
cnt = 0;

line = fgetl(fid_2D);
while ~all(line == -1) && cnt < cnt_max
    %fprintf([line, '\n'])
    
    if str_in_cell(line, timesteps)
        cnt_timestep = ind_str_in_cell(line, timesteps);
        data(cnt_timestep).n_horizon = sscanf(line, 'n_horizon: %d');
    elseif str_in_cell(line, maps)
        cnt_map = ind_str_in_cell(line, maps);
        data(cnt_timestep).map(cnt_map).map = line;
    elseif str_in_cell(line, solvers)
        cnt_solver = ind_str_in_cell(line, solvers);
        data(cnt_timestep).map(cnt_map).solver(cnt_solver).solver = line;
    else
        line_data = str2num(line);
        try
            data(cnt_timestep).map(cnt_map).solver(cnt_solver).trial_data = ...
                [data(cnt_timestep).map(cnt_map).solver(cnt_solver).trial_data; line_data];
        catch
            data(cnt_timestep).map(cnt_map).solver(cnt_solver).trial_data = line_data;
        end
    end
    
    if ~on_3D
        line = fgetl(fid_2D);
        if all(line == -1)
            on_3D = true;
        end
    end

    if on_3D
        line = fgetl(fid_3D);
    end

    cnt = cnt+1;
end

fclose(fid_2D);
fclose(fid_3D);

%% compute averages across trials

for ii = 1:length(data)
    for mm = 1:length(data(ii).map)
        for jj = 1:length(data(ii).map(mm).solver)
            data(ii).map(mm).solver(jj).avg_time = mean(data(ii).map(mm).solver(jj).trial_data(:,1));
            data(ii).map(mm).solver(jj).max_time = mean(data(ii).map(mm).solver(jj).trial_data(:,2));
            data(ii).map(mm).solver(jj).avg_nodes = mean(data(ii).map(mm).solver(jj).trial_data(:,3));
            data(ii).map(mm).solver(jj).max_nodes = mean(data(ii).map(mm).solver(jj).trial_data(:,4));
            data(ii).map(mm).solver(jj).qp_solve_time = mean(data(ii).map(mm).solver(jj).trial_data(:,5));
            data(ii).map(mm).solver(jj).avg_time_warmstart = mean(data(ii).map(mm).solver(jj).trial_data(:,6));
            data(ii).map(mm).solver(jj).max_time_warmstart = mean(data(ii).map(mm).solver(jj).trial_data(:,7));
            data(ii).map(mm).solver(jj).avg_nodes_warmstart = mean(data(ii).map(mm).solver(jj).trial_data(:,8));
            data(ii).map(mm).solver(jj).max_nodes_warmstart = mean(data(ii).map(mm).solver(jj).trial_data(:,9));
            data(ii).map(mm).solver(jj).qp_solve_time_warmstart = mean(data(ii).map(mm).solver(jj).trial_data(:,10));
        end
    end
end

%% plot avg and max solution times for 15 step case

% indexing for 15 steps
ind_15 = 3;

% get average and max times
avg_times = zeros(length(maps), length(solvers)); % init
max_times = zeros(length(maps), length(solvers)); % init
for ii = 1:length(maps)
    for jj = 1:length(solvers)
        if strcmp(solvers{jj}, 'ours_hz_warmstart')
            avg_times(ii,jj) = data(ind_15).map(ii).solver(jj).avg_time_warmstart;
            max_times(ii,jj) = data(ind_15).map(ii).solver(jj).max_time_warmstart;
        else
            avg_times(ii,jj) = data(ind_15).map(ii).solver(jj).avg_time;
            max_times(ii,jj) = data(ind_15).map(ii).solver(jj).max_time;
        end
    end
end

% bar plot options
num_decimals = 2;
do_60_clip = true;

% create figure
fig = figure();
til = tiledlayout(2,1);

% plot avg
max_val = 6;

nexttile
h = plotBar(maps_plot, avg_times, max_val, do_60_clip, false, face_colors, face_alphas, num_decimals);
title('Average Solution Time', 'FontSize', 14, 'Interpreter', 'latex')
ylim([0, 1.2*max_val])

% plot max
max_val = 20;

nexttile
h = plotBar(maps_plot, max_times, max_val, do_60_clip, true, face_colors, face_alphas, num_decimals);
title('Maximum Solution Time', 'FontSize', 14, 'Interpreter', 'latex')
ylim([0, 1.2*max_val])

% tile logic
ylabel(til, 'Solution Time [sec]', 'FontSize', 14, 'Interpreter', 'latex')
xlabel(til, 'Obstacle Map', 'FontSize', 14, 'interpreter', 'latex')

leg = legend(h, solvers_plot, ...
    'Orientation', 'horizontal', 'FontSize', 12, 'interpreter', 'latex');
leg.Layout.Tile = 'north';

til.TileSpacing = 'compact';
%til.Padding = 'tight';
til.Padding = 'compact';
fig.Position = [50 50 1180 490];

%% plot average and max b+b iterations

% indexing for 15 steps
ind_15 = 3;

% get average and max iterations
avg_nodes = zeros(length(maps), 3); % init
max_nodes = zeros(length(maps), 3); % init
for ii = 1:length(maps)
    for jj = 1:3
        % force cases where solver was not able to converge in 60 seconds
        % to NaN
        if data(ind_15).map(ii).solver(jj).avg_time > 60 || data(ind_15).map(ii).solver(jj).max_time > 60
            avg_nodes(ii,jj) = NaN;
            max_nodes(ii,jj) = NaN;
        elseif strcmp(solvers{jj}, 'ours_hz_warmstart')
            avg_nodes(ii,jj) = data(ind_15).map(ii).solver(jj).avg_nodes_warmstart;
            max_nodes(ii,jj) = data(ind_15).map(ii).solver(jj).max_nodes_warmstart;
        else
            avg_nodes(ii,jj) = data(ind_15).map(ii).solver(jj).avg_nodes;
            max_nodes(ii,jj) = data(ind_15).map(ii).solver(jj).max_nodes;
        end
    end
end

% create figure
fig = figure();
til = tiledlayout(2,1);

% bar plot options
num_decimals = 0;
do_60_clip = false;

% plot avg
max_val = 1.1*max(max(avg_nodes));

nexttile
h = plotBar(maps_plot, avg_nodes, max_val, do_60_clip, false, face_colors, face_alphas, num_decimals);
title('Average Branch and Bound Iterations', 'FontSize', 14, 'interpreter', 'latex')
ylim([0, 1.2*max_val])

% plot max
max_val = 1.1*max(max(max_nodes));

nexttile
h = plotBar(maps_plot, max_nodes, max_val, do_60_clip, true, face_colors, face_alphas, num_decimals);
title('Maximum Branch and Bound Iterations', 'FontSize', 14, 'interpreter', 'latex')
ylim([0, 1.2*max_val])

% tile logic
ylabel(til, 'Iterations', 'FontSize', 14, 'interpreter', 'latex')
xlabel(til, 'Obstacle Map', 'FontSize', 14, 'interpreter', 'latex')

leg = legend(h, solvers_plot{1:3}, ...
    'Orientation', 'horizontal', 'FontSize', 12, 'interpreter', 'latex');
leg.Layout.Tile = 'north';

til.TileSpacing = 'compact';
%til.Padding = 'tight';
til.Padding = 'compact';
fig.Position = [100 100 700 490];

%% plot average QP solution times

% indexing for 15 steps
ind_15 = 3;

% qp solve time
avg_qp = zeros(length(maps), 2); % init
for ii = 1:length(maps)
    for jj = 1:3
        if strcmp(solvers{jj}, 'ours_hz_warmstart')
            avg_qp(ii,jj) = data(ind_15).map(ii).solver(jj).qp_solve_time_warmstart*1000;
        else
            avg_qp(ii,jj) = data(ind_15).map(ii).solver(jj).qp_solve_time*1000;
        end
    end
end

% bar plot options
num_decimals = 0;
do_60_clip = false;

% max value
max_val = 200;

% create figure
fig = figure();
til = tiledlayout(1,1);

% plot avg
nexttile
h = plotBar(maps_plot, avg_qp, max_val, do_60_clip, true, face_colors, face_alphas, num_decimals);
title('Average QP Solution Times', 'FontSize', 14, 'interpreter', 'latex')
ylim([0, 1.2*max_val])

% tile logic
ylabel(til, 'Solution Times [ms]', 'FontSize', 14, 'interpreter', 'latex')
xlabel(til, 'Obstacle Map', 'FontSize', 14, 'interpreter', 'latex')

leg = legend(h, solvers_plot{1:3}, ...
    'Orientation', 'horizontal', 'FontSize', 12, 'interpreter', 'latex');
leg.Layout.Tile = 'north';

til.TileSpacing = 'compact';
%til.Padding = 'tight';
til.Padding = 'compact';
fig.Position = [100 100 700 420];

%% plot solution times vs number of iterations within a single map

for ii = 1:length(maps)

    % get average and max times
    avg_times = zeros(length(timesteps), length(solvers)); % init
    max_times = zeros(length(timesteps), length(solvers)); % init
    for ind = 1:length(timesteps)
        for jj = 1:length(solvers)
            if strcmp(solvers{jj}, 'ours_hz_warmstart')
                avg_times(ind,jj) = data(ind).map(ii).solver(jj).avg_time_warmstart;
                max_times(ind,jj) = data(ind).map(ii).solver(jj).max_time_warmstart;
            else
                avg_times(ind,jj) = data(ind).map(ii).solver(jj).avg_time;
                max_times(ind,jj) = data(ind).map(ii).solver(jj).max_time;
            end
        end
    end

    % bar plot options
    num_decimals = 2;
    do_60_clip = true;

    % create figure
    fig = figure();
    til = tiledlayout(2,1);

    % plot avg
    max_val = min([5, max(max(avg_times))]);

    nexttile
    h = plotBar(timesteps_num_str, avg_times, max_val, do_60_clip, false, face_colors, face_alphas, num_decimals);
    title('Average Solution Time', 'FontSize', 14, 'Interpreter', 'latex')
    ylim([0, 1.2*max_val])

    % plot max
    max_val = min([10, max(max(max_times))]);

    nexttile
    h = plotBar(timesteps_num_str, max_times, max_val, do_60_clip, true, face_colors, face_alphas, num_decimals);
    title('Maximum Solution Time', 'FontSize', 14, 'Interpreter', 'latex')
    ylim([0, 1.2*max_val])

    % tile logic
    ylabel(til, 'Solution Time [sec]', 'FontSize', 14, 'Interpreter', 'latex')
    xlabel(til, 'MPC Horizon Length $N$', 'FontSize', 14, 'interpreter', 'latex')

    leg = legend(h, solvers_plot, ...
        'Orientation', 'horizontal', 'FontSize', 12, 'interpreter', 'latex');
    leg.Layout.Tile = 'north';

    til.TileSpacing = 'compact';
    %til.Padding = 'tight';
    til.Padding = 'compact';
    fig.Position = [50 50 1180 490];

    % name window to keep track of maps
    fig.Name = ['map ', maps_plot{ii}];

end

%% read the no reach / no diag data

% baseline - copied over from other sim run for consistency
for cnt_timestep = 1:length(timesteps)
    data_bl(cnt_timestep).n_horizon = timesteps{cnt_timestep};
    
    cnt_map = ind_str_in_cell(line, maps);
    for cnt_map = 1:length(maps)
        data_bl(cnt_timestep).map(cnt_map).map = maps{cnt_map};

        % get solver index
        cnt_solver = ind_str_in_cell('ours_hz', solvers);

        % copy data over from other sim run
        data_bl(cnt_timestep).map(cnt_map).avg_time = data(cnt_timestep).map(cnt_map).solver(cnt_solver).avg_time;
        data_bl(cnt_timestep).map(cnt_map).max_time = data(cnt_timestep).map(cnt_map).solver(cnt_solver).max_time;
        data_bl(cnt_timestep).map(cnt_map).avg_nodes = data(cnt_timestep).map(cnt_map).solver(cnt_solver).avg_nodes;
        data_bl(cnt_timestep).map(cnt_map).max_nodes = data(cnt_timestep).map(cnt_map).solver(cnt_solver).max_nodes;
        data_bl(cnt_timestep).map(cnt_map).qp_solve_time = data(cnt_timestep).map(cnt_map).solver(cnt_solver).qp_solve_time;
    end
end

% no diag
fid_2D = fopen('solver_nodiag.txt', 'r');
fid_3D = fopen('solver_nodiag_3D.txt', 'r');
on_3D = false;
line = fgetl(fid_2D);
while ~all(line == -1)
    if str_in_cell(line, timesteps)
        cnt_timestep = ind_str_in_cell(line, timesteps);
        data_nd(cnt_timestep).n_horizon = sscanf(line, 'n_horizon: %d');
    elseif str_in_cell(line, maps)
        cnt_map = ind_str_in_cell(line, maps);
        data_nd(cnt_timestep).map(cnt_map).map = line;
    else
        line_data = str2num(line);
        try
            data_nd(cnt_timestep).map(cnt_map).trial_data = ...
                [data_nd(cnt_timestep).map(cnt_map).trial_data; line_data];
        catch
            data_nd(cnt_timestep).map(cnt_map).trial_data = line_data;
        end
    end

    if ~on_3D
        line = fgetl(fid_2D);
        if all(line == -1)
            on_3D = true;
        end
    end

    if on_3D
        line = fgetl(fid_3D);
    end
end
fclose(fid_2D);
fclose(fid_3D);

% no reach
fid_2D = fopen('solver_noreach.txt', 'r');
fid_3D = fopen('solver_noreach_3D.txt', 'r');
on_3D = false;
line = fgetl(fid_2D);
while ~all(line == -1)
    if str_in_cell(line, timesteps)
        cnt_timestep = ind_str_in_cell(line, timesteps);
        data_nr(cnt_timestep).n_horizon = sscanf(line, 'n_horizon: %d');
    elseif str_in_cell(line, maps)
        cnt_map = ind_str_in_cell(line, maps);
        data_nr(cnt_timestep).map(cnt_map).map = line;
    else
        line_data = str2num(line);
        try
            data_nr(cnt_timestep).map(cnt_map).trial_data = ...
                [data_nr(cnt_timestep).map(cnt_map).trial_data; line_data];
        catch
            data_nr(cnt_timestep).map(cnt_map).trial_data = line_data;
        end
    end

    if ~on_3D
        line = fgetl(fid_2D);
        if all(line == -1)
            on_3D = true;
        end
    end

    if on_3D
        line = fgetl(fid_3D);
    end
end
fclose(fid_2D);
fclose(fid_3D);

% no diag no reach
fid_2D = fopen('solver_nodiag_noreach.txt', 'r');
fid_3D = fopen('solver_nodiag_noreach_3D.txt', 'r');
on_3D = false;
line = fgetl(fid_2D);
while ~all(line == -1)
    if str_in_cell(line, timesteps)
        cnt_timestep = ind_str_in_cell(line, timesteps);
        data_ndnr(cnt_timestep).n_horizon = sscanf(line, 'n_horizon: %d');
    elseif str_in_cell(line, maps)
        cnt_map = ind_str_in_cell(line, maps);
        data_ndnr(cnt_timestep).map(cnt_map).map = line;
    else
        line_data = str2num(line);
        try
            data_ndnr(cnt_timestep).map(cnt_map).trial_data = ...
                [data_ndnr(cnt_timestep).map(cnt_map).trial_data; line_data];
        catch
            data_ndnr(cnt_timestep).map(cnt_map).trial_data = line_data;
        end
    end

    if ~on_3D
        line = fgetl(fid_2D);
        if all(line == -1)
            on_3D = true;
        end
    end

    if on_3D
        line = fgetl(fid_3D);
    end
end
fclose(fid_2D);
fclose(fid_3D);

% 1 thread
fid_2D = fopen('solver_1_thread.txt', 'r');
fid_3D = fopen('solver_1_thread_3D.txt', 'r');
on_3D = false;
line = fgetl(fid_2D);
while ~all(line == -1)
    if str_in_cell(line, timesteps)
        cnt_timestep = ind_str_in_cell(line, timesteps);
        data_1th(cnt_timestep).n_horizon = sscanf(line, 'n_horizon: %d');
    elseif str_in_cell(line, maps)
        cnt_map = ind_str_in_cell(line, maps);
        data_1th(cnt_timestep).map(cnt_map).map = line;
    else
        line_data = str2num(line);
        try
            data_1th(cnt_timestep).map(cnt_map).trial_data = ...
                [data_1th(cnt_timestep).map(cnt_map).trial_data; line_data];
        catch
            data_1th(cnt_timestep).map(cnt_map).trial_data = line_data;
        end
    end

    if ~on_3D
        line = fgetl(fid_2D);
        if all(line == -1)
            on_3D = true;
        end
    end

    if on_3D
        line = fgetl(fid_3D);
    end
end
fclose(fid_2D);
fclose(fid_3D);

% 8 thread
fid_2D = fopen('solver_8_thread.txt', 'r');
fid_3D = fopen('solver_8_thread_3D.txt', 'r');
on_3D = false;
line = fgetl(fid_2D);
while ~all(line == -1)
    if str_in_cell(line, timesteps)
        cnt_timestep = ind_str_in_cell(line, timesteps);
        data_8th(cnt_timestep).n_horizon = sscanf(line, 'n_horizon: %d');
    elseif str_in_cell(line, maps)
        cnt_map = ind_str_in_cell(line, maps);
        data_8th(cnt_timestep).map(cnt_map).map = line;
    else
        line_data = str2num(line);
        try
            data_8th(cnt_timestep).map(cnt_map).trial_data = ...
                [data_8th(cnt_timestep).map(cnt_map).trial_data; line_data];
        catch
            data_8th(cnt_timestep).map(cnt_map).trial_data = line_data;
        end
    end

    if ~on_3D
        line = fgetl(fid_2D);
        if all(line == -1)
            on_3D = true;
        end
    end

    if on_3D
        line = fgetl(fid_3D);
    end
end
fclose(fid_2D);
fclose(fid_3D);

% compute averages
for ii = 1:length(data_nr)
    for mm = 1:length(data_nr(ii).map)
        data_nr(ii).map(mm).avg_time = mean(data_nr(ii).map(mm).trial_data(:,1));
        data_nr(ii).map(mm).max_time = mean(data_nr(ii).map(mm).trial_data(:,2));
        data_nr(ii).map(mm).avg_nodes = mean(data_nr(ii).map(mm).trial_data(:,3));
        data_nr(ii).map(mm).max_nodes = mean(data_nr(ii).map(mm).trial_data(:,4));
        data_nr(ii).map(mm).qp_solve_time = mean(data_nr(ii).map(mm).trial_data(:,5));
    end
end

for ii = 1:length(data_nd)
    for mm = 1:length(data_nd(ii).map)
        data_nd(ii).map(mm).avg_time = mean(data_nd(ii).map(mm).trial_data(:,1));
        data_nd(ii).map(mm).max_time = mean(data_nd(ii).map(mm).trial_data(:,2));
        data_nd(ii).map(mm).avg_nodes = mean(data_nd(ii).map(mm).trial_data(:,3));
        data_nd(ii).map(mm).max_nodes = mean(data_nd(ii).map(mm).trial_data(:,4));
        data_nd(ii).map(mm).qp_solve_time = mean(data_nd(ii).map(mm).trial_data(:,5));
    end
end

for ii = 1:length(data_ndnr)
    for mm = 1:length(data_ndnr(ii).map)
        data_ndnr(ii).map(mm).avg_time = mean(data_ndnr(ii).map(mm).trial_data(:,1));
        data_ndnr(ii).map(mm).max_time = mean(data_ndnr(ii).map(mm).trial_data(:,2));
        data_ndnr(ii).map(mm).avg_nodes = mean(data_ndnr(ii).map(mm).trial_data(:,3));
        data_ndnr(ii).map(mm).max_nodes = mean(data_ndnr(ii).map(mm).trial_data(:,4));
        data_ndnr(ii).map(mm).qp_solve_time = mean(data_ndnr(ii).map(mm).trial_data(:,5));
    end
end

for ii = 1:length(data_1th)
    for mm = 1:length(data_1th(ii).map)
        data_1th(ii).map(mm).avg_time = mean(data_1th(ii).map(mm).trial_data(:,1));
        data_1th(ii).map(mm).max_time = mean(data_1th(ii).map(mm).trial_data(:,2));
        data_1th(ii).map(mm).avg_nodes = mean(data_1th(ii).map(mm).trial_data(:,3));
        data_1th(ii).map(mm).max_nodes = mean(data_1th(ii).map(mm).trial_data(:,4));
        data_1th(ii).map(mm).qp_solve_time = mean(data_1th(ii).map(mm).trial_data(:,5));
    end
end

for ii = 1:length(data_1th)
    for mm = 1:length(data_8th(ii).map)
        data_8th(ii).map(mm).avg_time = mean(data_8th(ii).map(mm).trial_data(:,1));
        data_8th(ii).map(mm).max_time = mean(data_8th(ii).map(mm).trial_data(:,2));
        data_8th(ii).map(mm).avg_nodes = mean(data_8th(ii).map(mm).trial_data(:,3));
        data_8th(ii).map(mm).max_nodes = mean(data_8th(ii).map(mm).trial_data(:,4));
        data_8th(ii).map(mm).qp_solve_time = mean(data_8th(ii).map(mm).trial_data(:,5));
    end
end

%% plot avg and max solution times for 15 step case - no reachability and no diagonal matrices

% indexing for 15 steps
ind_15 = 3;

% data
data_arr = {data_bl, data_nr, data_nd, data_ndnr};

% face colors and alphas
face_colors{1} = [0,0,1];
face_colors{2} = [255,179,0]/255;
face_colors{3} = [13,198,0]/255;
face_colors{4} = [212,33,99]/255;

face_alphas = {1,1,1,1};


% get average and max times
avg_times = zeros(length(maps), length(data_arr)); % init
max_times = zeros(length(maps), length(data_arr)); % init
for ii = 1:length(maps)
    for jj = 1:length(data_arr)
        avg_times(ii,jj) = data_arr{jj}(ind_15).map(ii).avg_time;
        max_times(ii,jj) = data_arr{jj}(ind_15).map(ii).max_time;
    end
end

% bar plot options
num_decimals = 2;
do_60_clip = true;

% create figure
fig = figure();
til = tiledlayout(2,1);

% plot avg
max_val = 1.5;

nexttile
h = plotBar(maps_plot, avg_times, max_val, do_60_clip, false, face_colors, face_alphas, num_decimals);
title('Average Solution Time', 'FontSize', 14, 'Interpreter', 'latex')
ylim([0, 1.2*max_val])

% plot max
max_val = 11;

nexttile
h = plotBar(maps_plot, max_times, max_val, do_60_clip, true, face_colors, face_alphas, num_decimals);
title('Maximum Solution Time', 'FontSize', 14, 'Interpreter', 'latex')
ylim([0, 1.2*max_val])

% tile logic
ylabel(til, 'Solution Time [sec]', 'FontSize', 14, 'Interpreter', 'latex')
xlabel(til, 'Obstacle Map', 'FontSize', 14, 'interpreter', 'latex')

leg = legend(h, {'Baseline', 'No Reachability', 'No Diagonal Matrix', 'No Reachability or Diagonal Matrix'}, ...
    'Orientation', 'horizontal', 'FontSize', 12, 'interpreter', 'latex');
leg.Layout.Tile = 'north';

til.TileSpacing = 'compact';
%til.Padding = 'tight';
til.Padding = 'compact';
fig.Position = [50 50 810 490];




%% embedded functions
function out = str_in_cell(str, cell_arr)

    out = false;
    for ii = 1:length(cell_arr)
        if strcmp(cell_arr{ii}, str)
            out = true;
            return;
        end
    end

end

function ind = ind_str_in_cell(str, cell_arr)
    
    ind = -1;
    for ii = 1:length(cell_arr)
        if strcmp(cell_arr{ii}, str)
            ind = ii;
            return;
        end
    end
    
end

function bp = plotBar(names, vals, max_val, do_60_clip, labels_on, face_colors, face_alphas, num_decimals)
    
    % clip according to max_val
    vals_clipped = vals;
    vals_clipped(isnan(vals_clipped)) = max_val;
    vals_clipped(vals_clipped > max_val) = max_val;
   
    % plot bars
    bp = bar(vals_clipped, 0.8);
    hold on;
    set(gca, 'XTick', 1:length(names))
    if labels_on
        set(gca, 'XTickLabel', names)
    else
        set(gca, 'XTickLabel', {[]})
    end
    set(gca, 'FontSize', 12)
    set(gca, 'GridAlpha', 0.05)
    set(gca , 'TickLabelInterpreter', 'latex')

    % color pattern
    for ii = 1:length(bp)
        bp(ii).FaceColor = face_colors{ii};
        bp(ii).FaceAlpha = face_alphas{ii};
    end

    % add number above bar
    bar_label_linespec = ['$%.', num2str(num_decimals), 'f$'];
    [n_maps,n_solvers] = size(vals);
    for ii = 1:n_maps
        for jj = 1:n_solvers

            % top of bar
            x_text = bp(jj).XEndPoints(ii);
            y_text = vals_clipped(ii,jj) + 0.025*max_val;
            star_on = false; % default
            if (vals(ii,jj) > 60 && do_60_clip) || isnan(vals(ii,jj))
                str = 'N/A';
                star_on = true;
            elseif vals(ii,jj) > max_val
                str = sprintf(bar_label_linespec, vals(ii,jj));
                star_on = true;
            else
                str = sprintf(bar_label_linespec, vals(ii,jj));
            end

            text(x_text, y_text, str, 'vert','bottom','horiz','center', 'FontSize', 12, 'interpreter', 'latex');

            % put star on bar if applicable
            if star_on
                % top of bar
                x_star = bp(jj).XEndPoints(ii);
                y_star = max_val;

                % star
                plot(x_star, y_star, '*k', 'MarkerSize', 10);
            end

        end
    end
end