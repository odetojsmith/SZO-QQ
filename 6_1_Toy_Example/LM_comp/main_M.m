clear all
M_set = [3 6 10 15 20];
L = 5;
i = 1;
[cost,timing,X_diff] = run_experiment(L,M_set(i));
COST_M = cost;
TIMING_M = timing;
X_DIFF_M = X_diff;

save('COST_M.mat','COST_M');
save('TIMING_M.mat','TIMING_M');
save('X_DIFF_M.mat','X_DIFF_M');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
M_set = [3 6 10 15 20];
L = 5;
i = 2;
[cost,timing,X_diff] = run_experiment(L,M_set(i));
load('COST_M.mat');
load('TIMING_M.mat');
load('X_DIFF_M.mat');
COST_M = [COST_M;cost];
TIMING_M = [TIMING_M;timing];
X_DIFF_M = [X_DIFF_M;X_diff];

save('COST_M.mat','COST_M');
save('TIMING_M.mat','TIMING_M');
save('X_DIFF_M.mat','X_DIFF_M');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
M_set = [3 6 10 15 20];
L = 5;
i = 3;
[cost,timing,X_diff] = run_experiment(L,M_set(i));
load('COST_M.mat');
load('TIMING_M.mat');
load('X_DIFF_M.mat');
COST_M = [COST_M;cost];
TIMING_M = [TIMING_M;timing];
X_DIFF_M = [X_DIFF_M;X_diff];

save('COST_M.mat','COST_M');
save('TIMING_M.mat','TIMING_M');
save('X_DIFF_M.mat','X_DIFF_M');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
M_set = [3 6 10 15 20];
L = 5;
i = 4;
[cost,timing,X_diff] = run_experiment(L,M_set(i));
load('COST_M.mat');
load('TIMING_M.mat');
load('X_DIFF_M.mat');
COST_M = [COST_M;cost];
TIMING_M = [TIMING_M;timing];
X_DIFF_M = [X_DIFF_M;X_diff];

save('COST_M.mat','COST_M');
save('TIMING_M.mat','TIMING_M');
save('X_DIFF_M.mat','X_DIFF_M');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
M_set = [3 6 10 15 20];
L = 5;
i = 5;
[cost,timing,X_diff] = run_experiment(L,M_set(i));
load('COST_M.mat');
load('TIMING_M.mat');
load('X_DIFF_M.mat');
COST_M = [COST_M;cost];
TIMING_M = [TIMING_M;timing];
X_DIFF_M = [X_DIFF_M;X_diff];

save('COST_M.mat','COST_M');
save('TIMING_M.mat','TIMING_M');
save('X_DIFF_M.mat','X_DIFF_M');



