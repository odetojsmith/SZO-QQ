clc; clear; close all;
%% Network specification
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
param.V_scale = 200;

% Load packages
% run('/Users/bguo/Library/Mobile Documents/com~apple~CloudDocs/Research Paper/zero-th order optimization sum/power system experiments/packages/matpower7.1/startup2.m');
addpath('/Users/ywang104/Documents/Master/Thesis/IFAC/pglib-opf-master')
addpath('/Users/ywang104/Documents/Master/Thesis/IFAC/pglib-opf-master/api')
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
[y, fi, result_next] = simulation(x_p, x_v,mpc,param);
result_current = result_next;
opf_f = runopf(result_current,mpopt).f;

dim = 2 * length(PV_bus)+1;
m = length(fi);     % Number of constraints
M0 = 0.5; % M_obj
Ms = 0.5; % M_con

d = dim;
M = M0;

param.dim = dim;
param.m = m;
param.M = M;
param.pvnumber = length(PV_bus);



%% Algorithm specification



reg= 0.0001;  % truncate alpha if alpha is too small; function: compute_gamma

n = 1; % n_sample for gradient estimate
nu = 0.001;  % step size for gradient estimate

sigma = 0.00; % std of noise
hat_sigma=sigma; % std of noise
eta = 0.0005; % termination condition on gradient |g|<eta;

delta = 0.1; % confidence level
T = 30000; % n_iteration
K = 1; % n_iteration outer loop
factor = 0.7; % decaying factor for eta at each outer loop
T_total = 0;

time_SGD = 0;

%% Run optimization
x = [x_p';x_v'];
x_progress = 1;
y_progress = 1;

rng('default')
s = rng;
rng(s)

epsilon = 5*1e-2;
th_convergence = 1e-8;
timing = 0;


[f_0, f_i] = func_realization(x,mpc,param);
f_opt  = 800.14;


rec_x = x;
rec_cost= f_0;
rec_max_con = max(f_i);


tic
for k=1:K

    xt = x;
    Tk = 0;
    xs = [];
    ys = [];
    fis = [];

    for t = 1:T
        [objective_value, constraints_values, objective_grad, constraints_grad, alphas, hat_sigma,s] =sample(xt,d,m,sigma,hat_sigma,delta,nu,n,mpc,param);
        step = dB_estimator(m,eta,alphas,constraints_grad,objective_grad,reg);
        step_norm = norm(step);
        gamma = compute_gamma(step,alphas,eta,delta,hat_sigma,constraints_grad,m,n,reg,M0,Ms);

        if step_norm < eta
            break;
        end
        xt = xt - gamma * step;
        if t == 1
            x_trajectory = xt;
            gamma_trajectory = gamma;
            [f_0, f_i] = func_realization(xt,mpc,param);

            errors_trajectory = f_0;
            constraints_trajectory = max(f_i);
            worst_constraint = max(f_i);
        else
            [f_0, f_i] = func_realization(xt,mpc,param);
            error = f_0;
            constraints = max(f_i);
            worst_constraint = max(f_i);

            x_trajectory = [x_trajectory xt];
            gamma_trajectory =[gamma_trajectory; gamma];
            errors_trajectory = [errors_trajectory, error];
            constraints_trajectory = [constraints_trajectory, constraints];
            worst_constraint = max(worst_constraint, worst_constraint);
        end
        xs = [xs xt];
        ys = [ys f_0];
        fis = [fis f_i];
        x_last = xt;
        fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", t,ys(end),max(fis(:,end)));

        time_SGD = [time_SGD; toc];
    end
    toc
    rec_cost = [rec_cost, errors_trajectory];
    rec_max_con = [rec_max_con constraints_trajectory];
    rec_x = [rec_x x_trajectory];
    T_total = T_total + Tk;
    x0 = x_last;
    eta = eta * factor;
end

rec_x_LB = rec_x;
rec_cost_LB = rec_cost;
rec_max_con_LB = rec_max_con;
rec_timing_LB = time_SGD';

save('rec_x_LB.mat','rec_x_LB');
save('rec_cost_LB.mat','rec_cost_LB');
save('rec_max_con_LB.mat','rec_max_con_LB')
save('rec_timing_LB.mat','rec_timing_LB');


function [objective_value, constraints_values, objective_grad, constraints_grad, alphas, hat_sigma,s ]=sample(x,d,m,sigma,hat_sigma,delta,nu,n,mpc,param)
% objective_value = f(x) + (sigma / n^0.5)*randn(1,1);
% constraints_values = h(x) + (sigma / n^0.5)*randn(m,1); 

[objective_value, constraints_values] = func_realization(x,mpc,param);

hat_sigma = sigma;

for j = 1:n
    s_unnormalized = randn(d, 1);
    s = s_unnormalized / norm(s_unnormalized);
    if j == 1
        [f_pert, h_pert] = func_realization(x + nu * s,mpc,param);
        objective_grad = (d*(f_pert + sigma*randn(1,1) - objective_value)/nu) * s ./ n;
        constraints_grad = ((((d *(h_pert + sigma*randn(m, 1) -constraints_values) ./ nu)* s')) ./ n);
    else
        objective_grad = objective_grad + (d *(f_pert+ sigma*randn(1,1) - objective_value)/nu) * s ./ n;
        constraints_grad = constraints_grad + ((((d *(h_pert + sigma*randn(m, 1) -constraints_values) ./ nu)* s')) ./ n);
    end
    alphas = - constraints_values -(log(1/ delta))^0.5 * sigma / n^0.5 * ones(m,1) - nu * ones(m,1);

end
end

function dB=dB_estimator(m,eta,alphas,constraints_grad,objective_grad,reg)

jacobian = constraints_grad;
df_e = objective_grad;
denominators = 1./max(ones(m,1)*reg, alphas);
dB = df_e + eta * jacobian'*denominators;

end

function gamma = compute_gamma(step,alphas,eta,delta,hat_sigma,constraints_grad,m,n,reg,M0,Ms)

step_norm = norm(step);
alphas_reg = alphas;

dhs = constraints_grad;
L_dirs = zeros(m,1);

for i = 1:m
    L_dirs(i) = abs(dot(dhs(i,:),step/step_norm))+sqrt(log(1/delta))*hat_sigma/sqrt(n);
    alphas_reg(i) = max(reg,alphas(i));
end
alphas = alphas_reg;
M2 = M0 + 2*eta*sum(Ms./alphas_reg) + 4 * eta * sum(L_dirs.^2 ./ alphas_reg.^2);
gamma = min(1 / step_norm * min(alphas ./ ( 2 * L_dirs +  alphas_reg.^0.5 * Ms^0.5)),1/M2 );

end





function [metric_value_obj, metric_value_cons] = func_realization(x,mpc,param)
[x_in,x_out] = split_state(x,param);
[y,fi,~] = simulation(x_in,x_out,mpc,param);
metric_value_obj = y;
metric_value_cons = fi';
end


function icon_pos = func_proxy(M,xi,l_max,delta_cons,metric_value_con_only)
delta_modi = delta_cons * xi * l_max+ 0.5*M * xi^2 * l_max^2; 
if max(metric_value_con_only + delta_modi)<0
    icon_pos = 0;
else
    icon_pos = 1;
end
end


function [y,fi,result] = simulation(x_in,x_out,mpc,param)
x_out = x_out/param.V_scale;
PG = 2;
VG = 6;
VM = 8;
PMAX = 9;
PV_bus = find(mpc.bus(:,2)==2);
REF_bus = 1;
valid_pg = find(mpc.gen(:,PMAX)~=0);
valid_pv = valid_pg(2:end);
mpopt = mpoption('verbose',0,'out.all',0,'opf.use_vg',0);
mpc.gen(:,VG) = x_out(1,:)';
mpc.bus([REF_bus',PV_bus'],VM)=x_out(1,:)';
mpc.gen(valid_pv,PG) = x_in(1,:)';
%     mpc.bus(PV_bus,BS) = x_in(1,2:end)';
result = runpf(mpc,mpopt);

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
