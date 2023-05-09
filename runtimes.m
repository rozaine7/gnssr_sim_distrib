%
% Collect statistics on the run times for various ensemble sizes
%
nicmat = [];
simrunmat = [];
comprunmat = [];
simcpumat = [];
compcpumat = [];

load runtimes100;
nicmat = [nicmat; nic];
simrunmat = [simrunmat; mean(simtime)];
comprunmat = [comprunmat; mean(comptime)];
simcpumat = [simcpumat; mean(simcpu)];
compcpumat = [compcpumat; mean(compcpu)];

load runtimes500;
nicmat = [nicmat; nic];
simrunmat = [simrunmat; mean(simtime)];
comprunmat = [comprunmat; mean(comptime)];
simcpumat = [simcpumat; mean(simcpu)];
compcpumat = [compcpumat; mean(compcpu)];

load runtimes1000;
nicmat = [nicmat; nic];
simrunmat = [simrunmat; mean(simtime)];
comprunmat = [comprunmat; mean(comptime)];
simcpumat = [simcpumat; mean(simcpu)];
compcpumat = [compcpumat; mean(compcpu)];

load runtimes5000;
nicmat = [nicmat; nic];
simrunmat = [simrunmat; mean(simtime)];
comprunmat = [comprunmat; mean(comptime)];
simcpumat = [simcpumat; mean(simcpu)];
compcpumat = [compcpumat; mean(compcpu)];

