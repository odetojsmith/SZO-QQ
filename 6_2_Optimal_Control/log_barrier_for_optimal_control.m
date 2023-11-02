clear all
%% Initialization
L=20;
M=20;

%%%%%%%%%%%%%% Computation of the initial point
mat_perf.Q = 5*eye(2); 
mat_perf.R = 2 * eye(2);
md.A = [1.1 1;-0.5 1.1];
md.B = eye(2);
md_p.A = [1.1 1;-0.5 1.1];
md_p.B = eye(2);
time.T = 6;
time.t = 1;
up_bound = 1.2;
x_0_oc = [1;1];
[K_0,~,~] = dlqr(md_p.A,md_p.B,mat_perf.Q,mat_perf.R,0);
K_0 = -K_0;
[cost,x_rec,u_rec] = sys_real_K(x_0_oc,time,md,K_0,mat_perf);
[cost,x_rec] = sys_real_u(x_0_oc,time,md,u_rec,mat_perf);
x_0_original = u_rec;

mat_perf.Q = 0.5*eye(2); 

%%%%%%%%%%%%%%%%%%%%% Log barrier for optimal control

eta = 0.00001;                 % Barrier parameter
T = 35000;                    % Number of iterations

TH = 1e-3;                 % Convergence condition

x0 = x_0_original';                % Starting point
% x0 = [0,-4.9];             % Second options, it is near a boundary 
y0 = obj_fun(x0);           % Initial objective
fi0 = fi_fun(x0,md);           % Initial constraint    
x_measurement = x0;         % x record
y_measurement = y0;         % y record
fi_measurement = fi0;       % fi record

gt_measurement = [];        % gradient record

d = size(x0,2);             % dimensions
m = size(fi0,2);            % number of constraints
timing_log = 0;
cost_log = obj_fun(x0);
x_extended = [obj_fun(x0)+0.000001;x0'];
x_log =  x_extended;
tic

%% Algorithm
for iter = 1:T
    x_current = x_measurement(end,:);
    y_current = y_measurement(end,:);
    fi_current = fi_measurement(end,:);
    vt = min(eta/(sqrt(d)*M),min(-fi_current)/max(L,sqrt(d)*M));

    X_current = repmat(x_current,d,1);
    VT = eye(d)*vt;
    X_current = X_current+VT;
    Y_current = [];
    FI_current = [];
    for i = 1:d
        Y_current = [Y_current; obj_fun(X_current(i,:))]; 
        FI_current = [FI_current; fi_fun(X_current(i,:),md)];
    end
    G0 = (Y_current - y_current)/vt;    % G0 = dx1
    GI = (FI_current - fi_current)/vt;  % GI = dxm

    gt = G0 + eta*sum(GI./(-fi_current),2);
    gt_measurement = [gt_measurement;gt];
    
    L2 = M + sum(2*eta*M./(-fi_current) + 4*eta*L^2./((-fi_current).^2));

    gama = min(min(-fi_current)/(2*L*norm(gt)),1/L2);

    
    
    
    x_next = x_current - gama*gt';
    y_next = obj_fun(x_next);
    fi_next = fi_fun(x_next,md);
    
    x_extended = [obj_fun(x_next)+0.000001;x_next'];
    x_log = [x_log x_extended];
    timing_log = [timing_log toc];
    cost_log = [cost_log obj_fun(x_next)];
    
    x_measurement = [x_measurement;x_next];
    y_measurement = [y_measurement;y_next];
    fi_measurement = [fi_measurement;fi_next];
    
    lambda = eta./-fi_next;

    %fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", iter,y_next,max(fi_next));
    
    
end


