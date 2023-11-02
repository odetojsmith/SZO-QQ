clear all 
L = 5;
M = 3;
[cost_log,timing_log,n_sample_log] = run_experiment_log(L,M);

save('cost_log.mat','cost_log')
save('timing_log.mat','timing_log')
save('n_sample_log.mat','n_sample_log')