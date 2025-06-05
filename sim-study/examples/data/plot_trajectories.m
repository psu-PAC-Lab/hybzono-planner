%% setup

clear all
close all

% flag: animate?
animate = false;

% map names and labels
maps = {'poly2', 'nonconvpoly1', 'zonopsu', 'zonocost1', 'zono3D'};
map_lbls = {'(a)', '(b)', '(c)', '(d)', '(e)'};

% start/end points
start_pts = {[-20,0], [-20,0], [-37.5,0], [0,2.5], [-15,0,0]};
end_pts = {[20,0], [20,0], [37.5,0], [25,2.5], [15,0,0]};

%% read text file

% 2D
fid = fopen('solver_comparison_traj.txt', 'r');

cnt_max = 1e4;
cnt = 0;

cnt_map = 0;

line = fgetl(fid);
while ~all(line == -1) && cnt < cnt_max
    
    if str_in_cell(line, maps)
        cnt_map = cnt_map+1;
        data2D(cnt_map).map = line;
    else
        line_data = str2num(line);
        try
            data2D(cnt_map).traj = ...
                [data2D(cnt_map).traj; line_data];
        catch
            data2D(cnt_map).traj = line_data;
        end
    end
    
    line = fgetl(fid);
    cnt = cnt+1;
end

fclose(fid);

% 3D
fid = fopen('solver_comparison_traj_3D.txt', 'r');

cnt = 0;
cnt_map = 0;

line = fgetl(fid);
while ~all(line == -1) && cnt < cnt_max
    
    if str_in_cell(line, maps)
        cnt_map = cnt_map+1;
        data3D(cnt_map).map = line;
    else
        line_data = str2num(line);
        try
            data3D(cnt_map).traj = ...
                [data3D(cnt_map).traj; line_data];
        catch
            data3D(cnt_map).traj = line_data;
        end
    end
    
    line = fgetl(fid);
    cnt = cnt+1;
end

% load motion plan
if animate
    
    % 2D
    fid = fopen('solver_comparison_motion_plan.txt', 'r');
    
    line = fgetl(fid);
    while ~all(line == -1)
        if str_in_cell(line, maps)
            cnt_map = find_in_cell(line, maps);
        elseif ~isempty(sscanf(line, 'k = %d'))
            k = sscanf(line, 'k = %d');
        else
            try
                data2D(cnt_map).motion_plan{k+1} = [data2D(cnt_map).motion_plan{k+1};
                                                    str2num(line)];
            catch
                data2D(cnt_map).motion_plan{k+1} = str2num(line);
            end
        end
        line = fgetl(fid);
    end

    fclose(fid);

    % 3D
    fid = fopen('solver_comparison_motion_plan_3D.txt', 'r');
    
    line = fgetl(fid);
    while ~all(line == -1)
        if str_in_cell(line, maps)
            cnt_map = 1; % lazy
        elseif ~isempty(sscanf(line, 'k = %d'))
            k = sscanf(line, 'k = %d');
        else
            try
                data3D(cnt_map).motion_plan{k+1} = [data3D(cnt_map).motion_plan{k+1};
                                                    str2num(line)];
            catch
                data3D(cnt_map).motion_plan{k+1} = str2num(line);
            end
        end
        line = fgetl(fid);
    end

    fclose(fid);

end

%% make 2D figure

n_maps_2D = length(data2D);

fig_2D = figure();
fig_2D.Position = [100, 100, 800, 450];

til_2D = tiledlayout(1, n_maps_2D);

for ii = 1:n_maps_2D

    % increment
    nexttile
    hold on;

    % get map
    map_file = [data2D(ii).map, '_vrep.txt'];
    [V_lim, V_obs, V_free, cost_vec] = read_map(map_file);

    % plot bounds
    V_lim_aug = [V_lim; V_lim(1,:)];
    plot(V_lim_aug(:,2), V_lim_aug(:,1), '-k');

    % plot obstacles
    for jj = 1:length(V_obs)
        patch(V_obs{jj}(:,2), V_obs{jj}(:,1), 'k')
    end

    % plot cell costs if applicable
    if strcmp(data2D(ii).map, 'zonocost1')

        % color map
        red_hsv = rgb2hsv([1,0,0]);
        cmap_hsv = repmat(red_hsv, 256, 1);
        cmap_hsv(:,2) = linspace(0,1,256)';
        cmap = hsv2rgb(cmap_hsv);

        min_cost = min(cost_vec);
        max_cost = max(cost_vec);

        % plot cost regions
        for jj = 1:length(cost_vec)
            frac_cost = (cost_vec(jj)-min_cost)/(max_cost-min_cost);
            int_cost = floor(frac_cost*256) + 1; % +1 for matlab indexing
            if (int_cost > 256) int_cost = 256; end
            p = patch(V_free{jj}(:,2), V_free{jj}(:,1), cmap(int_cost,:));
            p.EdgeAlpha = 0.2;
        end

        % colorbar
        colormap(cmap)
        colorbar
        clim([min_cost, max_cost])

    else % plot free space polytopes in white
        for jj = 1:length(V_free)
            p = patch(V_free{jj}(:,2), V_free{jj}(:,1), [1, 1, 1]);
            p.EdgeAlpha = 0.2;
        end
    end

    % plot start and end
    ind_map = find_in_cell(data2D(ii).map, maps);
    plot(start_pts{ind_map}(2), start_pts{ind_map}(1), 'Marker', 'square', 'MarkerFaceColor', [1,0,0], 'MarkerSize', 8, 'LineStyle', 'none', 'Color', [1,0,0])
    plot(end_pts{ind_map}(2), end_pts{ind_map}(1), 'Marker', 'pentagram', 'MarkerFaceColor', [0,1,0], 'MarkerSize', 8, 'LineStyle', 'none', 'Color', [0,1,0])

    % plot trajectory
    plot(data2D(ii).traj(:,2), data2D(ii).traj(:,1), '.b')

    % title
    title(['\textbf{', map_lbls{ind_map}, '}'], 'interpreter', 'latex', 'FontSize', 14);

    % plot limits
    ylim([start_pts{ind_map}(1), end_pts{ind_map}(1)])

end

% labels, etc
xlabel(til_2D, '$\eta$', 'interpreter', 'latex', 'fontsize', 14)
ylabel(til_2D, '$\zeta$', 'interpreter', 'latex', 'fontsize', 14)
til_2D.TileSpacing = 'compact';
til_2D.Padding = 'compact';

% make animation
if animate
    
    % white background
    set(gcf,'Color','w');

    % gif parameters
    gif_filename = 'sim_study_animated_trajectories_2D.gif';
    gif_delaytime = 0.3; % time between frames

    % number of time steps
    n_timesteps_map = zeros(1,length(data2D)); % init
    for ii = 1:length(data2D)
        n_timesteps_map(ii) = length(data2D(ii).motion_plan)-1;
    end
    n_timesteps = max(n_timesteps_map);

    % loop through timesteps and maps and animate
    first = true;
    for k = 1:n_timesteps

        h_cnt = 1;
        for ii = 1:length(data2D)

            % get timestep
            if k > n_timesteps_map(ii)
                k_map = n_timesteps_map(ii);
            else
                k_map = k;
            end

            % position plan and position
            pos_plan = data2D(ii).motion_plan{k_map}(:,[1,3]);
            pos = data2D(ii).traj(k_map,:);

            % plot plan and position
            nexttile(til_2D, ii)
            hp(h_cnt) = plot(pos_plan(:,2), pos_plan(:,1), 'or');
            h_cnt = h_cnt+1;
            hp(h_cnt) = plot(pos(2), pos(1), '.k', 'MarkerSize', 20);
            h_cnt = h_cnt+1;

        end

        % call gif
        if first
            gif(gif_filename, 'DelayTime', gif_delaytime, 'LoopCount', inf)
            first = false;
        else
            gif()
        end

        % delete graphics handle
        delete(hp)

    end

end


%% make 3D figure


n_maps_3D = length(data3D);

fig_3D = figure();
fig_3D.Position = [100, 100, 400, 350];

til_3D = tiledlayout(1, n_maps_3D);

for ii = 1:n_maps_3D

    % increment
    nexttile
    hold on;

    % get map
    map_file = [data3D(ii).map, '_vrep.txt'];
    [V_lim, V_obs, V_free] = read_map(map_file);

    % plot obstacles
    p_cnt = 1;
    for jj = 1:length(V_obs)
        P = Polyhedron(V_obs{jj});
        P.minHRep();
        facs = P.getFacet();
        x_verts = []; y_verts = []; z_verts = []; % init
        for kk = 1:length(facs)
            dim_range = max(facs(kk).V, [], 1) - min(facs(kk).V, [], 1);
            dim_flat = find(dim_range == min(dim_range));
            V_flat = facs(kk).V(:,setdiff(1:3,dim_flat));
            V_flat = orderPtsClockwise(V_flat);
            V = zeros(length(V_flat), 3); % init
            V(:,setdiff(1:3,dim_flat)) = V_flat;
            V(:,dim_flat) = facs(kk).V(:,dim_flat);
            x_verts = [x_verts, V(:,1)];
            y_verts = [y_verts, V(:,2)];
            z_verts = [z_verts, V(:,3)];
        end

        patch('XData', y_verts, 'YData', x_verts, 'ZData', z_verts, 'FaceColor', [0.5, 0.5, 0.5], 'FaceLighting', 'gouraud', 'FaceAlpha', 1.0, 'EdgeColor', 'none')
    end

    % plot start and end
    ind_map = find_in_cell(data3D(ii).map, maps);
    plot3(start_pts{ind_map}(2), start_pts{ind_map}(1), start_pts{ind_map}(3), 'Marker', 'square', 'MarkerFaceColor', [1,0,0], 'MarkerSize', 8, 'LineStyle', 'none', 'Color', [1,0,0])
    plot3(end_pts{ind_map}(2), end_pts{ind_map}(1), end_pts{ind_map}(3), 'Marker', 'pentagram', 'MarkerFaceColor', [0,1,0], 'MarkerSize', 8, 'LineStyle', 'none', 'Color', [0,1,0])

    % plot trajectory
    plot3(data3D(ii).traj(:,2), data3D(ii).traj(:,1), data3D(ii).traj(:,3), '.b')

    % apply bounds as plotting limits
    xlim([min(V_lim(:,2)), max(V_lim(:,2))]);
    ylim([min(V_lim(:,1)), max(V_lim(:,1))]);
    zlim([min(V_lim(:,3)), max(V_lim(:,3))]);

    % view angle and lighting
    view(60,50)
    light('Position',-[1,1,1]);
    lighting gouraud;

    % title
    title(['\textbf{', map_lbls{ind_map}, '}'], 'interpreter', 'latex', 'FontSize', 14);

    % labels
    xlabel('$\eta$', 'interpreter', 'latex', 'fontsize', 14)
    ylabel('$\zeta$', 'interpreter', 'latex', 'fontsize', 14)
    zlabel('$\gamma$', 'interpreter', 'latex', 'fontsize', 14)
end

% labels, etc
til_3D.TileSpacing = 'compact';
til_3D.Padding = 'compact';

% make animation
if animate
    
    % white background
    set(gcf,'Color','w');

    % gif parameters
    gif_filename = 'sim_study_animated_trajectories_3D.gif';
    gif_delaytime = 0.3; % time between frames

    % number of time steps
    n_timesteps_map = zeros(1,length(data3D)); % init
    for ii = 1:length(data3D)
        n_timesteps_map(ii) = length(data3D(ii).motion_plan)-1;
    end
    n_timesteps = max(n_timesteps_map);

    % loop through timesteps and maps and animate
    first = true;
    for k = 1:n_timesteps

        h_cnt = 1;
        for ii = 1:length(data3D)

            % get timestep
            if k > n_timesteps_map(ii)
                k_map = n_timesteps_map(ii);
            else
                k_map = k;
            end

            % position plan and position
            pos_plan = data3D(ii).motion_plan{k_map}(:,[1,3,5]);
            pos = data3D(ii).traj(k_map,:);

            % plot plan and position
            nexttile(til_3D, ii)
            hp(h_cnt) = plot3(pos_plan(:,2), pos_plan(:,1), pos_plan(:,3), 'or');
            h_cnt = h_cnt+1;
            hp(h_cnt) = plot3(pos(2), pos(1), pos(3), '.k', 'MarkerSize', 20);
            h_cnt = h_cnt+1;

        end

        % call gif
        if first
            gif(gif_filename, 'DelayTime', gif_delaytime, 'LoopCount', inf)
            first = false;
        else
            gif()
        end

        % delete graphics handle
        delete(hp)

    end

end



    
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

function ind = find_in_cell(str, cell_arr)
    ind = -1;
    for ii = 1:length(cell_arr)
        if strcmp(cell_arr{ii}, str)
            ind = ii;
            return;
        end
    end
end