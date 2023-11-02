clear all
L_set = [5 8 10 20 40];
M = 3;


i = 1;
[cost,timing,X_diff] = run_experiment(L_set(i),M);
COST_L = cost;
TIMING_L = timing;
X_DIFF_L = X_diff;

save('COST_L.mat','COST_L');
save('TIMING_L.mat','TIMING_L');
save('X_DIFF_L.mat','X_DIFF_L');

clear all
L_set = [5 8 10 20 40];
M = 3;
i = 2;
[cost,timing,X_diff] = run_experiment(L_set(i),M);

load('COST_L.mat');
load('TIMING_L.mat');
load('X_DIFF_L.mat');
COST_L = [COST_L;cost];
TIMING_L = [TIMING_L;timing];
X_DIFF_L = [X_DIFF_L;X_diff];
save('COST_L.mat','COST_L');
save('TIMING_L.mat','TIMING_L');
save('X_DIFF_L.mat','X_DIFF_L');

clear all
L_set = [5 8 10 20 40];
M = 3;
i = 3;
[cost,timing,X_diff] = run_experiment(L_set(i),M);

load('COST_L.mat');
load('TIMING_L.mat');
load('X_DIFF_L.mat');
COST_L = [COST_L;cost];
TIMING_L = [TIMING_L;timing];
X_DIFF_L = [X_DIFF_L;X_diff];
save('COST_L.mat','COST_L');
save('TIMING_L.mat','TIMING_L');
save('X_DIFF_L.mat','X_DIFF_L');

clear all
L_set = [5 8 10 20 40];
M = 3;
i = 4;
[cost,timing,X_diff] = run_experiment(L_set(i),M);

load('COST_L.mat');
load('TIMING_L.mat');
load('X_DIFF_L.mat');
COST_L = [COST_L;cost];
TIMING_L = [TIMING_L;timing];
X_DIFF_L = [X_DIFF_L;X_diff];
save('COST_L.mat','COST_L');
save('TIMING_L.mat','TIMING_L');
save('X_DIFF_L.mat','X_DIFF_L');

clear all
L_set = [5 8 10 20 40];
M = 3;
i = 5;
[cost,timing,X_diff] = run_experiment(L_set(i),M);

load('COST_L.mat');
load('TIMING_L.mat');
load('X_DIFF_L.mat');
COST_L = [COST_L;cost];
TIMING_L = [TIMING_L;timing];
X_DIFF_L = [X_DIFF_L;X_diff];
save('COST_L.mat','COST_L');
save('TIMING_L.mat','TIMING_L');
save('X_DIFF_L.mat','X_DIFF_L');