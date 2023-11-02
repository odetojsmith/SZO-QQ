clear; clc; close all
profile on;


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

V_scale = 200;
param.V_scale = V_scale;
% Initialization
[total_PD, total_QD] = total_load(mpc);
valid_pg = find(mpc.gen(:,PMAX)~=0);
valid_pv = valid_pg(2:end);
mpc.branch(:,RATE_A) = 600;
mpc.gen(valid_pg,PG) = (total_PD-100)/size(valid_pv,1);
vm = mpc.bus([REF_bus',PV_bus'],VM);
pg = mpc.gen(valid_pv,PG);
x_p = pg';
x_v = vm' * V_scale;
[~,fi,~] = simulation(x_p,x_v,mpc,param);

dim = 2 * length(PV_bus)+1;
m = length(fi);     % Number of constraints


param.dim = dim;
param.m = m;
param.pvnumber = length(PV_bus);





%% Compute the optimizer
x = [x_p';x_v'];
xhat = x;
[metric_value_cons_pert]= func_realization(x,mpc,param);
f0 = metric_value_cons_pert(1);
fi = metric_value_cons_pert(2:end);
opf_f = runopf(mpc,mpopt).f;


M = 0.13;
th_convergence = 1e-8;
timing = 0;

t = 1;
dt = 1;

% log barrier gain
log_gain = 0.1*ones(m,1);
y = -(f0-log_gain'*log10(-fi));


% configure es frequency components
A = ones(dim,1)*0.2;  % stable value  0.2
A(1:5) = ones(5,1)*0.2; % stable value    0.2
omega = 1*flip([1 sqrt(3) sqrt(5) sqrt(7) sqrt(9) sqrt(11) sqrt(13) sqrt(15) ...
    sqrt(17) sqrt(19) sqrt(21)]');
phase = ones(dim,1)*0;
int_gain = 0.3;   % stable value 0.05

% high pass filter
butterorder=1;
butterfreq=ones(dim,1)*0.4;  % stable value  0.4
butterfreq(1:5) = ones(5,1)*0.01; % stable value 0.01



% rec_filter_data
y_before = zeros(dim,butterorder+1)+y;
y_after = zeros(dim,butterorder+1);


%
delta_x = zeros(dim,1);

rec_x = [];
rec_obj = [];
rec_max_con = [];
rec_xhat = [];

% reduce_gain = false;
% reduce_gain2 = false;
% reduceA = 0;
tic
for k = 1:7000

    % change ES gain according to constraints
    if max(fi)>=-0.005 && ~reduce_gain2
        reduce_gain2 = true;
        int_gain = int_gain*0.7;
    end
    if max(fi)>=-0.002 && ~reduce_gain2
        reduce_gain2 = true;
        int_gain = int_gain*0.1;
    end

    % compute ES for every d
    for d = 1:dim
        [b,a] = butter(butterorder,butterfreq(d),'high');

        t = (k-1)*dt;

        [metric_value_cons_pert]= func_realization(x,mpc,param);
        f0 = metric_value_cons_pert(1);
        fi = metric_value_cons_pert(2:end);

        y = -(f0-log_gain'*log10(-fi));

        % high pass filter using difference equation
        y_before(d,1:butterorder) = y_before(d,2:butterorder+1);
        y_after(d,1:butterorder) = y_after(d,2:butterorder+1);

        y_before(d,butterorder+1) = y;

        y_new = 0;

        y_new = y_new + b*flip(y_before(d,:))';
        y_new = y_new - a(2)*y_after(d,1);
        y_new = y_new/a(1);
        y_after(d,2) = y_new;


        xi = y_new*sin(omega(d)*t + phase(d));
        delta_x(d) = xi*int_gain*dt;

    end

    xhat = xhat-delta_x;

    x = xhat + A.*sin(omega*t + phase);

    [metric_value_cons_pert]= func_realization(x,mpc,param);
    f0 = metric_value_cons_pert(1);
    fi = metric_value_cons_pert(2:end);


    rec_x = [rec_x x];
    rec_xhat = [rec_xhat xhat];
    rec_obj = [rec_obj f0];
    rec_max_con = [rec_max_con max(fi)];


    timing = [timing toc];

    fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f | Time: %1.2f \n", k,f0,max(fi),timing(end));
    k = k+1;

end

%%
timing = timing(1:end-1);
rec_cost_ES = rec_obj;
rec_max_con_ES = rec_max_con;
rec_timing_ES = timing;


save('rec_cost_ES.mat','rec_cost_ES');
save('rec_max_con_ES.mat','rec_max_con_ES');
save('rec_timing_ES.mat','rec_timing_ES');





function metric_value_cons = func_realization(x,mpc,param)
[x_in,x_out] = split_state(x,param);
[y,fi,~] = simulation(x_in,x_out,mpc,param);
%     y = y-x(1); %
%     if y>0
%         1+1;
%     end
metric_value_cons = [y;fi'];
end

function x_progress = line_search(M,x_dir,grad_cons,metric_value_cons,x,mpc,param)
delta_cons = (x_dir' * grad_cons(:,2:end))';
delta_pos = delta_cons(delta_cons>=10^(-7));
metric_value_con_only = metric_value_cons(2:end);
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
[metric_value_cons_ls]= func_realization(x_ls,mpc,param);
y_ls = metric_value_cons_ls(1);
y = metric_value_cons(1);

c = 1e-4;
while y_ls > y - c * xi * l_max * x_dir' * grad_cons(:,1)
    xi = xi * 0.6;
    x_ls = x + xi*l_max*x_dir;
    [metric_value_cons_ls]= func_realization(x_ls,mpc,param);
    y_ls = metric_value_cons_ls(1);

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

