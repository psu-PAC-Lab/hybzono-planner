close all
clear all

%% V-rep to hybzono

% get V-rep for whole map
map_file = 'poly2_vrep.txt';
[V_lim, V_obs, V_free, cost_vec] = read_map(map_file);

% plot whole map
for ii = 1:length(V_free)
    P(ii) = Polyhedron(V_free{ii});
end
ind_nonzero = find(P.volume > 1e-3);
P = P(ind_nonzero);

% plot
figure
plot(P)
hold on;
for ii = 1:length(P)
    poly_num_str = sprintf('%d', ii);
    xc = P(ii).chebyCenter.x;
    text(xc(1), xc(2), poly_num_str)
end

% take a subset of the polytopes
poly_nums = [3,4,11];

% get H-rep polytopes
P_old = P;
clear P
for ii = 1:length(poly_nums)
    P(ii) = P_old(poly_nums(ii));
end
P.computeHRep();


% get hybzono

% get V matrix
V = [];
for i = 1:length(P)
    V = [V; P(i).V];
end

ndec = 3;
V = round(V,ndec);
V = unique(V,'rows');

% get M matrix
M = zeros(size(V,1),length(P));
for i = 1:length(P)
    for j = 1:size(P(i).V,1)
        M(ismember(V,round(P(i).V(j,:),ndec),'rows')==1,i)=1;
    end
end

Zh = vPoly2hybZono({V', M});

% get relaxation of union of H-rep polytopes
bigM = 4; % 4 is smallest integer big enough to cover set
n_regions = length(P);
A_h0 = []; % init
b_h0 = []; % init
M_s = []; % init, selector matrix with big-M
for ii = 1:n_regions
    A_h0 = [A_h0;
        P(ii).A];
    b_h0 = [b_h0;
        P(ii).b + bigM*ones(size(P(ii).b))];

    M_s_ii = zeros(size(P(ii).A, 1), n_regions);
    M_s_ii(:,ii) = bigM*ones(size(P(ii).A, 1), 1);
    M_s = [M_s;
        M_s_ii];
end
n_cons = size(A_h0, 1);

A_h = [A_h0, M_s;
       zeros(n_regions, 2), -eye(n_regions);
       zeros(n_regions, 2), eye(n_regions);
       zeros(1, 2), -ones(1, n_regions);
       zeros(1, 2), ones(1, n_regions)];
b_h = [b_h0;
       zeros(n_regions,1);
       ones(n_regions,1);
       -1;
       1];

P_rel = Polyhedron(A_h, b_h);

P_rel_proj = projection(P_rel, [1,2], 'mplp'); % project onto spatial dimensions

% get relaxation of hybzono
Zc = conZono([Zh.Gc, Zh.Gb], Zh.c, [Zh.Ac, Zh.Ab], Zh.b);

% plot everything together
figure
hold on;
plot(P_rel_proj, 'color', [174,26,41]/255, 'alpha', 0.8);
plot(Zc, [101,62,191]/255, 0.8);
plot(P, 'color', [193,134,208]/255, 'alpha', 0.8);

h(1) = patch(NaN, NaN, [193,134,208]/255);
h(2) = patch(NaN, NaN, [101,62,191]/255);
h(3) = patch(NaN, NaN, [174,26,41]/255);
legend(h, {'$\mathcal{Z}$', '$CR(\mathcal{Z})$, HZ', '$CR(\mathcal{Z})$, H-Rep Big-M'}, 'interpreter', 'latex');

%% occupancy grid

clear all

% load map
x_range = [-37.5, 37.5];
y_range = [-10, 10];

[~, P, Zh] = zonotopeTiling_PSU(x_range, y_range);

% get relaxation of union of H-rep polytopes
bigM = 75; % 4 is smallest integer big enough to cover set
n_regions = length(P);
A_h0 = []; % init
b_h0 = []; % init
M_s = []; % init, selector matrix with big-M
for ii = 1:n_regions
    A_h0 = [A_h0;
        P(ii).A];
    b_h0 = [b_h0;
        P(ii).b + bigM*ones(size(P(ii).b))];

    M_s_ii = zeros(size(P(ii).A, 1), n_regions);
    M_s_ii(:,ii) = bigM*ones(size(P(ii).A, 1), 1);
    M_s = [M_s;
        M_s_ii];
end
n_cons = size(A_h0, 1);

A_h = [A_h0, M_s;
       zeros(n_regions, 2), -eye(n_regions);
       zeros(n_regions, 2), eye(n_regions);
       zeros(1, 2), -ones(1, n_regions);
       zeros(1, 2), ones(1, n_regions)];
b_h = [b_h0;
       zeros(n_regions,1);
       ones(n_regions,1);
       -1;
       1];

P_rel = Polyhedron(A_h, b_h);

P_rel_proj = projection(P_rel, [1,2], 'mplp'); % project onto spatial dimensions

% get relaxation of hybzono
Zc = conZono([Zh.Gc, Zh.Gb], Zh.c, [Zh.Ac, Zh.Ab], Zh.b);

% plot everything together
figure
hold on;
plot(P_rel_proj, 'color', [174,26,41]/255, 'alpha', 0.8);
plot(Zc, [101,62,191]/255, 0.8);
plot(P, 'color', [193,134,208]/255, 'alpha', 0.8);
axis equal

h(1) = patch(NaN, NaN, [193,134,208]/255);
h(2) = patch(NaN, NaN, [101,62,191]/255);
h(3) = patch(NaN, NaN, [174,26,41]/255);
legend(h, {'$\mathcal{Z}$', '$CR(\mathcal{Z})$, HZ', '$CR(\mathcal{Z})$, H-Rep Big-M'}, 'interpreter', 'latex');


