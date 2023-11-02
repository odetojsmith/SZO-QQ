%%
clear; clc; close all
profile on;

M = 0.5;     % Smoothness
nu = 0.001;  % Fixed step length for gradient approximation
mu = 0.0001;  % Quadratic penalty coeff
xi = 0.02;   % Terminal conditions
rho = 0.1;

global nu;
nu = 0.001;

infeasible = 0;

param.V_scale = 200;
rec_largest_con = [];
rec_cost = [];

% Constants for matpower
BASE_KV = 10;
BS = 6;
RATE_A = 6;
VMIN = 13;
VMAX = 12;
PMAX = 9;
PMIN = 10;
QMIN = 5;
QMAX = 4;
PG = 2;
VG = 6;
VM = 8;


% Load packages
run('../packages/matpower7.1/startup2.m');
addpath('../packages/pglib-opf-master')
addpath('../packages/pglib-opf-master/api')
mpopt = mpoption('verbose',0,'out.all',0,'opf.use_vg',0);
mpc = loadcase('pglib_opf_case30_as.m');

% Network Specifications
PV_bus = find(mpc.bus(:,2)==2);
REF_bus = 1;
mpc.bus(:,VMAX) = 1.1;
mpc.bus(:,VMIN) = 0.9;
mpc.gen(:,QMIN) = -300;
mpc.gen(:,QMAX) = 300;
mpc.gen(1:2,PMAX) = [300;300];




% Initialization
[total_PD, total_QD] = total_load(mpc);
valid_pg = find(mpc.gen(:,PMAX)~=0);
valid_pv = valid_pg(2:end);
mpc.branch(:,RATE_A) = 600;
mpc.gen(valid_pg,PG) = (total_PD-100)/size(valid_pv,1);
vm = mpc.bus([REF_bus',PV_bus'],VM);
pg = mpc.gen(valid_pv,PG);
x_p = pg';
x_v = vm' * param.V_scale;
param.mpc = mpc;

[y, fi, result_next] = simulation(x_p,x_v,param);
result_current = result_next;
opf_f = runopf(result_current,mpopt).f;

dim = 2 * length(PV_bus)+1;
m = length(fi);     % Number of constraints


param.dim = dim;
param.m = m;
param.M = M;
param.pvnumber = length(PV_bus);




%% Compute the optimizer
x = [x_p';x_v'];
x_progress = 1;
y_progress = 1;

rng('default')
s = rng;
rng(s)


epsilon = 5*1e-2;
th_convergence = 1e-9;
timing = 0;

x_rec = x;
rec_cost = y;
rec_largest_con = max(fi);



x_var = sdpvar(dim,1);
ops = sdpsettings('verbose',0,'debug',1);


icon_first = 0;
icon_term = 1;
t=1;

tic
while (norm(y_progress)>th_convergence) 
    % Derive the gradients and function values
    data  = data_update(x,param);
    % Use the data to form a sub linear programming problem
   [x_dir,epsilon] = sub_LP(data,epsilon,param);
    % Implement the line search
    x_progress = line_search(M,x_dir,data,x,param);
    % Move forward
    x = x+x_progress;


    % Data recording
    x_sim = x';
    timing = [timing toc];
    [y,fi,~] = simulation(x_sim(1:param.pvnumber),x_sim(param.pvnumber+1:end),param);
    x_rec = [x_rec x];
    rec_cost = [rec_cost y];
    rec_largest_con = [rec_largest_con max(fi)];

    if icon_first == 1
        y_progress = rec_cost(end)-rec_cost(end-1);
    end
    icon_first = 1;
    fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f | Time: %1.2f \n", t,y,max(fi),timing(end));
    t = t+1;

end
toc
%%

rec_cost_LP = rec_cost;
rec_max_con_LP = rec_largest_con;
rec_timing_LP = timing;

save('rec_cost_LP.mat','rec_cost_LP');
save('rec_max_con_LP.mat','rec_max_con_LP');
save('rec_timing_LP.mat','rec_timing_LP');




function metric_value_funcs = func_realization(x,param)
[x_in,x_out] = split_state(x,param);
[y,fi,~] = simulation(x_in,x_out,param);

metric_value_funcs = [y;fi'];
end

function x_progress = line_search(M,x_dir,data,x,param)
    delta_cons = (x_dir' * data.grad_funcs(:,2:end))';
    delta_pos = delta_cons(delta_cons>=10^(-7));
    metric_value_con_only = data.metric_value_funcs(2:end);
    metric_value_pos = abs(metric_value_con_only(delta_cons>=10^(-7)));
    l_max = min(metric_value_pos./delta_pos);
    xi_max = 1;
    xi_min = 0;
    for i = 1:20
        xi = (xi_max+xi_min)/2;
        icon_pos = func_proxy(M,xi,l_max,delta_cons,metric_value_con_only);
        if icon_pos == 0
            xi_min = xi;
        else
            xi_max = xi;
        end
    end

    % line search with func evaluation
    x_ls = x + xi*l_max*x_dir;
    [metric_value_funcs_ls]= func_realization(x_ls,param);
    y_ls = metric_value_funcs_ls(1);
    y = data.metric_value_funcs(1);
    
    c = 1e-4;
    while y_ls > y - c * xi * l_max * x_dir' * data.grad_funcs(:,1)
        xi = xi * 0.6;
        x_ls = x + xi*l_max*x_dir;
        [metric_value_funcs_ls]= func_realization(x_ls,param);
        y_ls = metric_value_funcs_ls(1);
    
    end
    x_progress = 0.9 * xi * l_max * x_dir;
end

function icon_pos = func_proxy(M,xi,l_max,delta_cons,metric_value_con_only)
    delta_modi = delta_cons * xi * l_max+ 0.5*M * xi^2 * l_max^2; % add 0.5*
    if max(metric_value_con_only + delta_modi)<0
        icon_pos = 0;
    else
        icon_pos = 1;
    end
end


function [y,fi,result] = simulation(x_in,x_out,param)
    x_out = x_out/param.V_scale;
    PG = 2;
    VG = 6;
    VM = 8;
    PMAX = 9;
    PV_bus = find(param.mpc.bus(:,2)==2);
    REF_bus = 1;
    valid_pg = find(param.mpc.gen(:,PMAX)~=0);
    valid_pv = valid_pg(2:end);
    mpopt = mpoption('verbose',0,'out.all',0,'opf.use_vg',0);
    param.mpc.gen(:,VG) = x_out(1,:)';
    param.mpc.bus([REF_bus',PV_bus'],VM)=x_out(1,:)';
    param.mpc.gen(valid_pv,PG) = x_in(1,:)';
    %     mpc.bus(PV_bus,BS) = x_in(1,2:end)';
    result = runpf(param.mpc,mpopt);
    
    loss = sum(totcost(result.gencost,result.gen(:,PG)));
    
    y = loss;
    fi = fi_fun(result);
end


function fi = fi_fun(mpc)
    RATE_A = 6;
    VM = 8;
    VMIN = 13;
    VMAX = 12;
    PMAX = 9;
    PMIN = 10;
    QMIN = 5;
    QMAX = 4;
    QG = 3;
    PG = 2;
    rate_a = mpc.branch(:,RATE_A);
    pf = vecnorm(mpc.branch(:,14:15),2,2);
    pt = vecnorm(mpc.branch(:,16:17),2,2);
    vm = mpc.bus(:,VM);
    f1 = (mpc.bus(:,VMIN) - vm)';
    f2 = (vm - mpc.bus(:,VMAX))';
    f3 = (pf - rate_a)';
    f4 = (pt - rate_a)';
    % f5 = (mpc.gen(1:2,PMIN) - mpc.gen(1:2,PG))';
    % f6 = (mpc.gen(1:2,PG) - mpc.gen(1:2,PMAX))';
    % f7 = (mpc.gen(:,QMIN) - mpc.gen(:,QG))';
    % f8 = (mpc.gen(:,QG) - mpc.gen(:,QMAX))';
    fi = [f1,f2,f3,f4];
    %fi = [f1,f2,f3,f4,f5,f6,f7,f8];
end

function [x_in,x_out] = split_state(x,param)
    x_in = x(1:param.pvnumber);
    x_out = x(param.pvnumber+1:end);
    x_in = x_in';
    x_out = x_out';
end

function [x_dir,epsilon] = sub_LP(data,epsilon,param)
    
    while true

        f = data.grad_funcs(:,1)';
        A = [];
        for i = 1 : param.m
            if data.metric_value_funcs(i+1) > -epsilon
                grad_funcs_nm = data.grad_funcs(:,i+1)/norm(data.grad_funcs(:,i+1));  % add normalization
                A = [A; grad_funcs_nm'];
            end
        end
        b = -1*epsilon*ones(size(A,1),1);
        lb = -1*ones(size(f,2),1);
        ub = 1*ones(size(f,2),1);

        options = optimoptions('linprog','Display','none');

        [x_dir,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,options);

        if exitflag ~=1
            disp('something went wrong');
        end
        if fval>=0
            disp('not a descent direction');
        end

        if x_dir' *  data.grad_funcs(:,1)> -epsilon
            epsilon = epsilon/10;
        else
            break;
        end
    end
end

function data  = data_update(x,param)
    global nu;
    data.metric_value_funcs = func_realization(x,param);
    grad_funcs = [];
    for j = 1 : param.dim
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        [metric_value_funcs_pert]= func_realization(x_copy,param);
        grad_funcs_j = (metric_value_funcs_pert-data.metric_value_funcs)/nu;
        grad_funcs = [grad_funcs;grad_funcs_j'];
    end
    data.grad_funcs = grad_funcs;
end

