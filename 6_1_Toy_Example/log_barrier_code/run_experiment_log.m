function [cost_log,timing_log,n_sample_log] = run_experiment_log(L,M)
%% Parameters

eta = 0.001;                 % Barrier parameter
T = 60000;                    % Number of iterations

TH = 1e-6;                 % Convergence condition

x0 = [0.9,0.9];                % Starting point
% x0 = [0,-4.9];             % Second options, it is near a boundary 
y0 = obj_fun(x0);           % Initial objective
fi0 = fi_fun(x0);           % Initial constraint    
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
count = 0;
n_sample_log = 1;

%% Algorithm
for iter = 1:T
    count = count + 1;
    x_current = x_measurement(end,:);
    y_current = y_measurement(end,:);
    fi_current = fi_measurement(end,:);
    vt = min(eta/(sqrt(d)*M),min(-fi_current)/max(L,sqrt(d)*M));

    X_current = repmat(x_current,d,1);
    VT = eye(d)*vt;
    X_current = X_current+VT;
    Y_current = [];
    FI_current = [];
    t_1 = toc;
    for i = 1:d
        Y_current = [Y_current; obj_fun(X_current(i,:))]; 
        FI_current = [FI_current; fi_fun(X_current(i,:))];
    end
    t_2 = toc;
    G0 = (Y_current - y_current)/vt;    % G0 = dx1
    GI = (FI_current - fi_current)/vt;  % GI = dxm

    gt = G0 + eta*sum(GI./(-fi_current),2);
    gt_measurement = [gt_measurement;gt];
    
    L2 = M + sum(2*eta*M./(-fi_current) + 4*eta*L^2./((-fi_current).^2));

    gama = min(min(-fi_current)/(2*L*norm(gt)),1/L2);

    
    
    
    x_next = x_current - gama*gt';
    y_next = obj_fun(x_next);
    fi_next = fi_fun(x_next);
    
    x_extended = [obj_fun(x_next)+0.000001;x_next'];
    x_log = [x_log x_extended];
    timing_log = [timing_log toc-(t_2-t_1)];
    cost_log = [cost_log obj_fun(x_next)];
    n_sample_log = [n_sample_log 3*count+1];
    
    x_measurement = [x_measurement;x_next];
    y_measurement = [y_measurement;y_next];
    fi_measurement = [fi_measurement;fi_next];
    
    lambda = eta./-fi_next;

    %fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", iter,y_next,max(fi_next));
    
    
end





