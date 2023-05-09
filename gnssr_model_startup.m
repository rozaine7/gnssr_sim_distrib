rootpath = '/Users/jgarriso/My Drive/jgarriso/projects/gnssr/gnssr_sim_distrib/';
path('/Users/jgarriso/jgarriso/gnss_basic1.4', path)
path(strcat(rootpath,'prepro1.3'), path)
path(strcat(rootpath,'pdfret1.1'), path)
path(strcat(rootpath,'gnssr_theory3.0'), path)
path(strcat(rootpath,'postpro1.0'), path)
path(strcat(rootpath,'compdata1.2'), path)
path(strcat(rootpath,'simulator0.1'), path)
%path('/home/bfrc/a/jgarriso/public-web/share/m_map', path)

%
% Uncomment to compile the mex files
%
mex -compatibleArrayDims ./gnssr_theory3.0/dp_dtheta_dopp_sum_mex.c ...
    ./gnssr_theory3.0/dp_dtheta_dopp_sumf.c ./gnssr_theory3.0/gram_charlierf.c ...
    -outdir ./gnssr_theory3.0/

%beep off