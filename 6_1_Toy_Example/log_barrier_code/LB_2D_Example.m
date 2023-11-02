clear; clc; close all
%% Parameters
L = 5;                      % Lipschitz
M = 3;                      % Smoothness
eta = 0.001;                 % Barrier parameter
T = 100000;                    % Number of iterations

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
        FI_current = [FI_current; fi_fun(X_current(i,:))];
    end
    G0 = (Y_current - y_current)/vt;    % G0 = dx1
    GI = (FI_current - fi_current)/vt;  % GI = dxm

    gt = G0 + eta*sum(GI./(-fi_current),2);
    gt_measurement = [gt_measurement;gt];
    
    L2 = M + sum(2*eta*M./(-fi_current) + 4*eta*L^2./((-fi_current).^2));

    gama = min(min(-fi_current)/(2*L*norm(gt)),1/L2);

    
    
    
    x_next = x_current - gama*gt'
    y_next = obj_fun(x_next);
    fi_next = fi_fun(x_next);
    
    timing_log = [timing_log toc];
    cost_log = [cost_log obj_fun(x_next)];
    
    x_measurement = [x_measurement;x_next];
    y_measurement = [y_measurement;y_next];
    fi_measurement = [fi_measurement;fi_next];
    
    lambda = eta./-fi_next;

    %fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", iter,y_next,max(fi_next));
    
    
    if abs((y_measurement(end,1)-y_measurement(end-1,1)))<=TH
        disp('converged');
        break
    end
end



function y = obj_fun(x)
    y = x(2)+0.1*x(1)^2;
end

function fi = fi_fun(x)
    f1 = -x(1);
    f2 = x(1)^2-x(2);
    f3 = x(2)-1;

    fi = [f1,f2,f3];
end

function []=plot_figure(x_measurement,y_measurement,fi_measurement)
    figure(1)
    plot(y_measurement(:,1),'b-'); hold on
    plot(-5*ones(size(y_measurement(:,1))),'r--'); hold on
    title('Objective');
    xlabel('$k$','Interpreter','latex')
    ylabel('$f^0(x)$','Interpreter','latex')
    legend('objective')

    figure(2)
    plot(fi_measurement(:,1),'b--'); hold on
    plot(fi_measurement(:,2),'b-.'); hold on
    plot(zeros(size(fi_measurement(:,1))),'r--'); hold on
    title('Constraints')
    legend('$x_1 \leq 2.7$','$x_2 \geq -5$','Interpreter','latex','Location','northwest')
    xlabel('$k$','Interpreter','latex')
    ylabel('$f^1(x)$','Interpreter','latex')

    figure(3)
    plot(x_measurement(:,1),x_measurement(:,2),'b-','MarkerFaceColor','b'); hold on
    
    line([2.7,2.7],[-5,1],'Color','r','Linewidth',1);
    line([-4,2.7],[-5,-5],'Color','r','Linewidth',1);
    scatter(x_measurement(1,1),x_measurement(1,2),80,'ro','filled','MarkerFaceColor','r')
    scatter(2.7,0.5,100,'gp','filled','MarkerFaceColor','g')

    title('Optimization trajectory')
    xlabel('$x_1$','Interpreter','latex')
    ylabel('$x_2$','Interpreter','latex')
    legend('0-LBM','Boundary 1','Boundary 2','Start','Optimum','Location','northwest')
    
    xlim([-4,3])
    ylim([-6,1])
    
end