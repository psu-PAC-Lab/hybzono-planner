function [V_lim, V_obs, V_free, cost_vec] = read_map(filename)

% read data
fid = fopen(filename,'r');

cnt_max = 1e4;
cnt = 0;
cnt_poly = 0;
obs_on = false;
free_on = false;
lim_on = false;
cost_on = false;
cost_vec = []; % init



line = fgetl(fid);
while ~all(line == -1) && cnt < cnt_max
    if strcmp(line, 'Obstacles')
        cnt_poly = 1;
        obs_on = true;
        free_on = false;
        lim_on = false;
        cost_on = false;
    elseif strcmp(line, 'Free Space')
        cnt_poly = 1;
        obs_on = false;
        free_on = true;
        lim_on = false;
        cost_on = false;
    elseif strcmp(line, 'Limits')
        obs_on = false;
        free_on = false;
        lim_on = true;
        cost_on = false;
    elseif strcmp(line, 'Costs')
        obs_on = false;
        free_on = false;
        lim_on = false;
        cost_on = true;
    elseif strcmp(line, '---')
        cnt_poly = cnt_poly + 1;
    else
        if obs_on
            try
                V_obs{cnt_poly} = [V_obs{cnt_poly}; str2num(line)];
            catch
                V_obs{cnt_poly} = str2num(line);
            end
        elseif free_on
            try
                V_free{cnt_poly} = [V_free{cnt_poly}; str2num(line)];
            catch
                V_free{cnt_poly} = str2num(line);
            end
        elseif lim_on
            try
                V_lim = [V_lim; str2num(line)];
            catch
                V_lim = str2num(line);
            end
        elseif cost_on
            try
                cost_vec = [cost_vec; str2num(line)];
            catch
                cost_vec = str2num(line);
            end
        end
    end

    line = fgetl(fid);
    cnt = cnt+1;
end

fclose(fid);

end