clear; clc; close all
profile on;

nu = 0.001;  % Fixed step length for gradient approximation
mu = 0.0001;    % Quadratic penalty coeff
xi = 0.002;  % Terminal conditions
% x_0 = [0.9;0.9];               % Initial point
param.V_scale = 50;
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

% run('/Users/bguo/Library/Mobile Documents/com~apple~CloudDocs/Research Paper/zero-th order optimization sum/power system experiments/packages/matpower7.1/startup2.m');
% addpath('/Users/bguo/Library/Mobile Documents/com~apple~CloudDocs/Research Paper/zero-th order optimization sum/power system experiments/packages/pglib-opf-master')
% addpath('/Users/bguo/Library/Mobile Documents/com~apple~CloudDocs/Research Paper/zero-th order optimization sum/power system experiments/packages/pglib-opf-master/api')
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
dim = 2 * length(PV_bus)+1;



% Initialization
[total_PD, total_QD] = total_load(mpc);
valid_pg = find(mpc.gen(:,PMAX)~=0);
valid_pv = valid_pg(2:end);
mpc.branch(:,RATE_A) = 600;
mpc.gen(valid_pg,PG) = (total_PD-100)/size(valid_pv,1);
vm = mpc.bus([REF_bus',PV_bus'],VM);
pg = mpc.gen(valid_pv,PG);
x0_in = pg';
x_out = vm' * param.V_scale;
[y0_in, fi0_in, result_next] = simulation(x0_in, x_out,mpc,param);
result_current = result_next;
opf_f = runopf(result_current,mpopt).f;
m = length(fi0_in);     % Number of constraints
M = 0.13;                  % Smoothness

param.dim = dim;
param.m = m;
param.M = M;
param.pvnumber = length(PV_bus);

% 


%% Compute the optimizer
x_0 = [y0_in+0.1;x0_in';x_out']; 
x_rec = x_0; 
x = x_0;
x_old = x_0+1;                  
x_diff = x-x_old;
i = 1;
while max(abs(x_diff(:))) >= xi
    grad_cons = [];
    CON = [];
    x_old = x;
    x_progress = sdpvar(dim+1,1);
    metric_value_cons = func_realization(x,mpc,param);
    for j = 1 : dim+1
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        metric_value_cons_pert = func_realization(x_copy,mpc,param);
        grad_cons_j = (metric_value_cons_pert-metric_value_cons)/nu;
        grad_cons = [grad_cons;grad_cons_j'];
    end
    if i == 1
        x_var = sdpvar(dim+1,1);
        grad_var = sdpvar(dim+1,m+1);
        metric_var = sdpvar(m+1,1);
        obj = [1 zeros(1,dim)] * x_progress + mu * (x_progress-x_var)' * (x_progress-x_var);

        for j = 1 : m+1
            CON_j = metric_var(j) + grad_var(:,j)' * ((x_progress-x_var)) + 2*M * (x_progress-x_var)'* (x_progress-x_var)<=0;
            CON = [CON, CON_j];
        end

        P = optimizer(CON,obj,[],{x_var,grad_var,metric_var},x_progress);
        x = P(x,grad_cons,metric_value_cons);
        break;
    end
    x = P(x,grad_cons,metric_value_cons);
    x_diff = x-x_old;
    i = i+1;
end

%% Start the iteration
x_0 = [y0_in+0.1;x0_in';x_out']; 
x_rec = x_0; 
x = x_0;
x_old = x_0+1;                  
x_diff = x-x_old;
i = 1;
timing = 0;
cost = f_0(x(2:end));

while max(abs(x_diff(:))) >= xi
    grad_cons = [];

    CON = [];
    x_old = x;
    x_progress = sdpvar(dim+1,1);
    metric_value_cons = func_realization(x,mpc,param);
    for j = 1 : dim+1
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        [metric_value_cons_pert]= func_realization(x_copy,mpc,param);
        grad_cons_j = (metric_value_cons_pert-metric_value_cons)/nu;
        grad_cons = [grad_cons;grad_cons_j'];
    end
    
   
    tic
    x = P(x,grad_cons,metric_value_cons);
    x_rec = [x_rec x];
    x_sim = x';
    timing = [timing timing(end)+toc];
    [y,fi,~] = simulation(x_sim(2:param.pvnumber+1),x_sim(param.pvnumber+2:end),mpc,param);
    rec_cost = [rec_cost y];
    rec_largest_con = [rec_largest_con max(fi)];
    x_diff = x-x_old;
    i = i+1;
    fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f | Time: %1.2f \n", i,y,max(fi),timing(end));
end

x = x_rec;
rec_cost_QQ = rec_cost;
rec_max_con_QQ = rec_largest_con;
rec_timing_QQ = timing(1:end-1);

opti = 800.14 * ones(size(rec_cost));
save('rec_cost_QQ.mat','rec_cost_QQ');
save('rec_max_con_QQ.mat','rec_max_con_QQ');
save('rec_timing_QQ.mat','rec_timing_QQ');






function y=f_0(x)
%%%%%%% Here enter the evaluation of the objective function
    y = x(2)+0.1*x(1)^2;
end

% function y=f_c(x,k)
%     if k == 1
%         y = x(1)-2;
%     end
%     if k == 2
%         y = -x(1);
%     end
%     if k == 3
%         y = x(2)-2;
%     end
%     if k == 4
%         y = -x(2);
%     end
% end

function y=f_c(x,k)
%%%%%%%%% Here, enter the evaluation of the constraint functions
    if k == 1
        y = -x(1);
    end
    if k == 2
        y = x(2)-1;
    end
    if k == 3
        y = x(1)^2-x(2);
    end
end

function metric_value_cons = func_realization(x,mpc,param)
    x_pure = x(2:end);
    [x_in,x_out] = split_state(x_pure,param);
    [y,fi,~] = simulation(x_in,x_out,mpc,param);
    y = y-x(1);
    if y>0
        1+1;
    end
    metric_value_cons = [y;fi'];
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

