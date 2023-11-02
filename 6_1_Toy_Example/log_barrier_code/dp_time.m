ETA = [1 0.9 0.8 0.7 0.6 0.5 0.35 0.3 0.2 0.17 0.16];
C = [];
m = 3;
M = 3;
L = 2;

for i = 1:length(ETA)
    C = [C 4*(m*M*ETA(i)+4*L^2*m)/(ETA(i))^(3)];
end

load('timing_log.mat')
C = ceil(C)
time_termi_rec = timing_log(C)
ETA_log = ETA;
save('ETA_log.mat','ETA_log');
save('time_termi_rec.mat','time_termi_rec');