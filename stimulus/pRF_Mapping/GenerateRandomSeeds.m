% generate random seeds 
% setting the seeds to random values is of course equivalent to not setting. But, this way I can figure out if
% I'm always using the same values (eg because the computer is always restarted at the same time before
% scanning...)

clc
nruns = 10; % how many scans we are running

% format bank
% disp(ceil(rand(nruns,1)*10^4))

for i = 1:nruns
    fprintf('%i\n',ceil(rand(1)*10^4))
end