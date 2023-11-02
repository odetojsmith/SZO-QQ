clear; clc; close all

% addpath('PATH: pglib-opf-master')
% addpath('PATH: matpower')

% Two loop optimization:
% Inner loop: optimize via active power generation (PG) 
% Outer loop: optimize via voltage magnitude of generator (VM)
% Once inner loop converges, go one step outer loop, until global convergence
%% Define LB parameters
eta = 0.01;
L_in = 0.5;
M_in = 1;
L_out = 10;
M_out = 50;  
T = 30000;

TH_in = 1e-2; % converge condition inner loop
TH = 1e-4;    % converge condition global
%% Define power flow model
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
mpopt = mpoption('verbose',0,'out.all',0,'opf.use_vg',0);
% mpc = loadcase('case9');
% mpc = loadcase('pglib_opf_case3_lmbd.m');
% mpc = loadcase('pglib_opf_case14_ieee.m');
mpc = loadcase('pglib_opf_case30_as.m');
% mpc = loadcase('pglib_opf_case118_ieee.m');
% mpc = loadcase('pglib_opf_case200_activ.m');

[total_PD, total_QD] = total_load(mpc);
total_PD = sum(total_PD);
[REF_bus, ~, PQ_bus] = bustypes(mpc.bus,mpc.gen);
PV_bus = find(mpc.bus(:,2)==2);
% Change box constraint limits
mpc.bus(:,VMAX) = 1.1;
mpc.bus(:,VMIN) = 0.9;
mpc.gen(:,QMIN) = -300;
mpc.gen(:,QMAX) = 300;
mpc.gen(1:2,PMAX) = [300;300];

valid_pg = find(mpc.gen(:,PMAX)~=0);    % Select buses whose PMAX not equals 0

valid_pv = valid_pg(2:end);             % Exclude slack bus
mpc.branch(:,RATE_A) = 600;

mpc.gen(valid_pg,PG) = (total_PD-100)/size(valid_pv,1); % Evenly assign initial PG

% vm = mpc.bus([REF_bus',PV_bus'],VM);
vm = mpc.gen(:,VG');
pg = mpc.gen(valid_pv,PG);

x_in = pg'; 
x_out = vm';
[y0_in, fi0_in, result_next] = simulation(x_in, x_out,mpc);
result_current = result_next;

x_measurement_in = x_in;
y_measurement_in = y0_in;
fi_measurement_in = fi0_in;


x_measurement_out = [];
y_measurement_out = [];
fi_measurement_out = [];

d_in = size(valid_pv,1);
d_out = size(PV_bus,1)+1;
m = size(fi0_in,2);

opf_f = runopf(result_current,mpopt).f;

%% Main loop
t = 0;
for iter = 1:T
    
    while true
        x_current_in = x_measurement_in(end,:);
        y_current_in = y_measurement_in(end,:);
        fi_current_in = fi_measurement_in(end,:);
        vt_in = min(eta/(sqrt(d_in)*M_in), min(-fi_current_in)/max(L_in,m*sqrt(d_in)*M_in)); %
    
        X_current_in = repmat(x_current_in,d_in,1);
        VT_in = eye(d_in)*vt_in;
        X_current_in = X_current_in+VT_in;
        Y_current_in = [];
        FI_current_in = [];
        for i = 1:d_in
            [y_e_in, fi_e_in,~] = simulation(X_current_in(i,:),x_out,result_current);
            Y_current_in = [Y_current_in; y_e_in]; 
            FI_current_in = [FI_current_in; fi_e_in];
        end
        G0_in = (Y_current_in - y_current_in)/vt_in;
        GI_in = (FI_current_in - fi_current_in)/vt_in;
    
        gt_in = G0_in + eta*sum(GI_in./(-fi_current_in),2); %;
        L2_in = M_in + sum(2*eta*M_in./(-fi_current_in) + 4*eta*L_in^2./((-fi_current_in).^2));
    
        gama_in = min(min(-fi_current_in)/(2*L_in*norm(gt_in)),1/L2_in);
    
        x_next_in = x_current_in - gama_in*gt_in';
        [y_next_in, fi_next_in, result_next] = simulation(x_next_in,x_out,result_current);
    
        x_measurement_in = [x_measurement_in;x_next_in];
        y_measurement_in = [y_measurement_in;y_next_in];
        fi_measurement_in = [fi_measurement_in; fi_next_in];
        result_current = result_next;
    
        lambda = eta./-fi_next_in;

        figure(1)
        plot(y_measurement_in,'b-','Linewidth',1);
        line([0,size(y_measurement_in,1)],[opf_f,opf_f],'Color','r','Linewidth',1,'LineStyle','--');
        legend('Generation cost','Benchmark');
        title('Objective inner loop');

        t = t+1;
        fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", t,y_next_in,max(fi_next_in));

        if abs((y_measurement_in(end,1)-y_measurement_in(end-1,1)))<=TH_in
            break;
        end
    end

if abs((y_measurement_in(end,1)-y_measurement_in(end-1,1)))<=TH
    break;
end
x_current_out = result_current.gen(:,VG)';
y_current_out = y_measurement_in(end,1);
fi_current_out = fi_measurement_in(end,:);

x_measurement_out = [x_measurement_out; x_current_out];
y_measurement_out = [y_measurement_out; y_current_out];
fi_measurement_out = [fi_measurement_out; fi_current_out];

vt_out = min(eta/(sqrt(d_out)*M_out), min(-fi_current_out)/max(L_out,m*sqrt(d_out)*M_out));
X_current_out = repmat(x_current_out,d_out,1);
VT_out = eye(d_out)*vt_out;
X_current_out = X_current_out+VT_out;
Y_current_out = [];
FI_current_out = [];
for i = 1:d_out
    [y_e_out, fi_e_out,~] = simulation(x_next_in,X_current_out(i,:),result_current);
    Y_current_out = [Y_current_out; y_e_out]; 
    FI_current_out = [FI_current_out; fi_e_out];
end
G0_out = (Y_current_out - y_current_out)/vt_out;
GI_out = (FI_current_out - fi_current_out)/vt_out;

gt_out = G0_out+ eta*sum(GI_out./(-fi_current_out),2);
L2_out = M_out + sum(2*eta*M_out./(-fi_current_out) + 4*eta*L_out^2./((-fi_current_out).^2));

gama_out = min(min(-fi_current_out)/(2*L_out*norm(gt_out)),1/L2_out);

x_next_out = x_current_out - gama_out*gt_out';
x_out = x_next_out;

figure(2)
plot(y_measurement_out,'b-',LineWidth=1); hold on
line([0,size(y_measurement_out,1)],[opf_f,opf_f],'Color','r','Linewidth',1,'LineStyle','--');
legend('Generation cost','Benchmark');
title('Objective outer loop');

t = t+1;
fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", t,y_next_in,max(fi_next_in));

end



function [y,fi,result] = simulation(x_in,x_out,mpc)
    PG = 2;
    VG = 6;
    VM = 8;
    PMAX = 9;
    [REF_bus, PV_bus, PQ_bus] = bustypes(mpc.bus,mpc.gen);
    PV_bus = find(mpc.bus(:,2)==2);
    valid_pg = find(mpc.gen(:,PMAX)~=0);
    valid_pv = valid_pg(2:end);
    mpopt = mpoption('verbose',0,'out.all',0,'opf.use_vg',0);
    mpc.gen(:,VG) = x_out(1,:)';
    mpc.bus([REF_bus',PV_bus'],VM)=x_out(1,:)';
    mpc.gen(valid_pv,PG) = x_in(1,:)';

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
    f1 = (pf - rate_a)';                         % branch power limits
    f2 = (pt - rate_a)';                         % branch power limits
    f3 = (mpc.bus(:,VMIN) - vm)';                % bus voltage magnitude limits
    f4 = (vm - mpc.bus(:,VMAX))';                % bus voltage magnitude limits
    f5 = (mpc.gen(1:2,PMIN) - mpc.gen(1:2,PG))'; % PG limits
    f6 = (mpc.gen(1:2,PG) - mpc.gen(1:2,PMAX))'; % PG limits
    f7 = (mpc.gen(:,QMIN) - mpc.gen(:,QG))';     % QG limits
    f8 = (mpc.gen(:,QG) - mpc.gen(:,QMAX))';     % QG limits
    fi = [f1,f2,f3,f4,f5,f6,f7,f8];
end
