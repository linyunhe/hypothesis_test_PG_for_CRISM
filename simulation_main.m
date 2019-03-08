% Constant-scalar simulation running
% N is the number of runs you'd like to do
% num_iter is the number of iterations you'd like to run
scalar = 10^4;
N = 10;
num_iter = 2500;
% Get correct path
cur_folder = pwd;
save_path = fullfile(cur_folder,'example_data/simulations_results');
sb_path = fullfile(cur_folder,'example_data/sb_and_wa/cdr490947778566_sb0000000l_3_envi.img');
wa_path = fullfile(cur_folder,'example_data/sb_and_wa/cdr490947778566_wa0000000l_3_envi.img');
ddr_file = fullfile(cur_folder,'example_data/ground_truth_ddr.img');

algorithm1_path = fullfile(cur_folder,'Algorithm1');
algorithm2_path = fullfile(cur_folder,'Algorithm2');


truth = read_envi_data(fullfile(cur_folder,'example_data/ground_truth_mu.bsq'));
for k = 1:N
    display(k)
% Generate scaled poisson distributed noise to the groundtruth    
    data = poissrnd(truth.*scalar)./scalar;
    iof = zeros(181,640,238);
    iof(1:170,32:631,:) = data;
    write_envi_data(iof,fullfile(save_path,sprintf('Simulations%i',k)),'bsq');
% Start to run algorithm1 and algorithm2
    cd(algorithm1_path); 

    ssa_file = fullfile(save_path,sprintf('Simulations%i.bsq',k));
    output_prefix = fullfile(save_path,sprintf('Simulations%i_PRUN',k));
    test5v5_compact(...
        'crism_iof_filename', ssa_file,...
        'ddr_iof_filename', ddr_file,...
        'sb_filename', sb_path,...
        'wa_filename', wa_path,...
        'num_iter', num_iter,...
        'switchover_at_iter', round(num_iter/5),...
        'switchover_flag',1,...
        'observation_base', output_prefix...        
    );

    cd(algorithm2_path)

    output_prefix = fullfile(save_path,sprintf('Simulations%i_GRUN',k));

    test5v5_compact_MSEE(...
        'crism_iof_filename', ssa_file,...
        'ddr_iof_filename', ddr_file,...
        'sb_filename', sb_path,...
        'wa_filename', wa_path,...
        'num_iter', num_iter,...
        'switchover_at_iter', round(num_iter/5),...
        'switchover_flag',1,...
        'observation_base', output_prefix...
        );

    cd(cur_folder)


end


% Compute the mean relative error of each runs
truth = read_envi_data(fullfile(cur_folder,'example_data/ground_truth_c.bsq'));
err_p = zeros(261,924,238,N);
err_g = zeros(261,924,238,N);
for k = 1:N
    est_p = read_envi_data(fullfile(save_path,sprintf('Simulations%i_PRUN_c_iter_%i.bsq',k,num_iter)));
    est_g = read_envi_data(fullfile(save_path,sprintf('Simulations%i_GRUN_c_iter_%i.bsq',k,num_iter)));
    err_p(:,:,:,k) = (est_p - truth)./truth;
    err_g(:,:,:,k) = (est_g - truth)./truth;
end
mean_err_p = mean(err_p,4);
mean_err_g = mean(err_g,4);
write_envi_data(mean_err_p,fullfile(save_path,sprintf('Simulations%iruns_mean_err_PRUN',N)),'bsq');
write_envi_data(mean_err_g,fullfile(save_path,sprintf('Simulations%iruns_mean_err_GRUN',N)),'bsq');
% mean relative error of PRUN
% mean(mean_err_p(~isnan(mean_err_p)))
% std relative error of PRUN
% sqrt(var(mean_err_p(~isnan(mean_err_p))))
% mean relative error of GRUN
% mean(mean_err_g(~isnan(mean_err_g)))
% std relative error of GRUN
% sqrt(var(mean_err_g(~isnan(mean_err_g))))