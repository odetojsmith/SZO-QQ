clear all
L = 5;
M_set = [3 6 10 15 20];
COST_M_log = [];
TIMING_M_log = [];
N_SAMPLE_M_log = [];
for i = 1: length(M_set)
    [cost,timing,n_sample] = run_experiment_log(L,M_set(i));
    COST_M_log = [COST_M_log;cost];
    TIMING_M_log = [TIMING_M_log;timing];
    N_SAMPLE_M_log = [N_SAMPLE_M_log; n_sample];
end

L_set = [5 8 10 20 40];
M = 3;
COST_L_log = [];
TIMING_L_log = [];
N_SAMPLE_L_log = [];

for i = 1: length(L_set)
    [cost,timing,n_sample] = run_experiment_log(L_set(i),M);
    COST_L_log = [COST_L_log;cost];
    TIMING_L_log = [TIMING_L_log;timing];
    N_SAMPLE_L_log = [N_SAMPLE_L_log; n_sample];
end

save('COST_M_log.mat','COST_M_log');
save('COST_L_log.mat','COST_L_log');
save('TIMING_M_log.mat','TIMING_M_log');
save('TIMING_L_log.mat','TIMING_L_log');



