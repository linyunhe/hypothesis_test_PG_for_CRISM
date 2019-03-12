function res = PG_model_selection(data_folder)
    cur_folder = pwd;
    algorithm1_path = fullfile(cur_folder,'Algorithm1');
    algorithm2_path = fullfile(cur_folder,'Algorithm2');
    params = load_options_file(fullfile(data_folder,'PG_running_info.txt'));
    res_save_file = fullfile(params.save_folder,'PG_running_results.txt');
    output_prefix_P = 'Poisson_assumped';
    output_prefix_G = 'Gaussian_assumped';
    num_iter_P = 30;
    num_iter_G = 100;
    use_all_rows = strcmp(params.rows_use_all,'true');
    cd(algorithm1_path); 
    test5v5_compact(...
        'crism_iof_filename', params.ssa_filename,...
        'ddr_iof_filename', params.ddr_filename,...
        'sb_filename', params.sb_filename,...
        'wa_filename', params.wa_filename,...
        'rows_use_all', use_all_rows,...
        'rows_subset_bounds', [str2num(params.row_min),str2num(params.row_max)],...
        'num_iter', num_iter_P,...
        'save_spectrogram_flag', 0,...
        'switchover_flag',0,...
        'bands_bounds',[str2num(params.band_min),str2num(params.band_max)],...
        'spatial_penalty_flag',0,...
        'spectral_penalty_flag',0,...
        'output_path_start_in',params.save_folder,...
        'cols_use_all',true,...
        'save_input_flag',true,...
        'observation_base', output_prefix_P...        
    );
    
   
    
    cd(algorithm2_path)
    test5v5_compact_MSEE(...
        'crism_iof_filename', params.ssa_filename,...
        'ddr_iof_filename', params.ddr_filename,...
        'sb_filename', params.sb_filename,...
        'wa_filename', params.wa_filename,...
        'rows_use_all', use_all_rows,...
        'rows_subset_bounds', [str2num(params.row_min),str2num(params.row_max)],...
        'num_iter', num_iter_G,...
        'switchover_flag',0,...
        'bands_bounds',[str2num(params.band_min),str2num(params.band_max)],...
        'spatial_penalty_flag',0,...
        'spectral_penalty_flag',0,...
        'output_path_start_in',params.save_folder,...
        'cols_use_all',true,...
        'observation_base', output_prefix_G...        
    );

    cd(cur_folder)
    ssa = read_envi_data(fullfile(params.save_folder,'Poisson_assumped_ssa.bsq'));
    mu_P = read_envi_data(fullfile(params.save_folder,sprintf('Poisson_assumped_mu_iter_%i.bsq',num_iter_P)));
    mu_G = read_envi_data(fullfile(params.save_folder,sprintf('Gaussian_assumped_mu_iter_%i.bsq',num_iter_G)));
    
    G2 =  2*( ssa(:).*log(ssa(:)./mu_P(:)) + ( mu_P(:) - ssa(:)));
    [d_kl_p,scalar_p] = Calculate_best_scalar(G2,10^4,10^5);
    
    varnorm = mean((ssa(:) - mu_G(:)).^2);
    G2_norm =  (ssa(:)-mu_G(:)).^2./(varnorm);
    pvalue = 1 - chi2cdf(G2_norm(:),1);
    [N,~] = hist(pvalue,100);
    pdf = N/sum(N);
    d_kl_g = sum(pdf.*log(pdf)) + log(100);
    
    pvalue_p = 1-chi2cdf(G2*scalar_p,1);
    p = [0:0.001:1];
    cdp = p;
    cdn = p;
%
    for i = 1:length(p)
        tmpp = (pvalue_p < p(i));
        tmp = (pvalue < p(i));
        cdp(i) = sum(tmpp)/length(pvalue_p);
        cdn(i) = sum(tmp)/length(pvalue);
    end
    
    figure;
    plot(p,p);
    hold on;
    plot(p,cdp,'o');
    hold on;
    plot(p,cdn,'x')
    title('Cumulative distribution of P-values');
    xlabel('Type I error \alpha');
    ylabel('Cumulative distribution of rejected measures');
    legend('Y=X','Scaled Poisson','Gaussian');
    savefig(fullfile(params.save_folder,'PG_running_fig.fig'));
    
    display('Writing results out...')
    fileID = fopen(res_save_file,'w');
    fprintf(fileID,'Poisson and Gaussian model selection results:\r\n');
    if d_kl_g > d_kl_p
        fprintf(fileID,'Scaled Poisson is selected!\r\n');
    else
        fprintf(fileID,'Gaussian is selected!\r\n');
    end
    fprintf(fileID,'Input Information: %s\r\n',fullfile(data_folder,'PG_running_info.txt'));
    fprintf(fileID,'Scaled Poisson:\n Scalar:%f\n d_kl:%f\r\n',scalar_p,d_kl_p);
    fprintf(fileID,'Gaussian:\n d_kl:%f\r\n',d_kl_g);    
    fclose(fileID);
    display('Writing results out...done!')
    res = true;
end


function params = load_options_file(opts_filename)
% read options from the given file. Expect exactly 1 option saved per line.
% Accepts literally any strings
fid = fopen(opts_filename);
fprintf('Loading options from %s...\n',opts_filename);
next_line = fgetl(fid);
while ~(isnumeric(next_line) && next_line == -1) % file has lines left
    result = regexp(next_line,'(\w+):\s*(.*)$','tokens');
    if ~isempty(result) % there is a match in this line
        param_name = result{1}{1};
        value = result{1}{2};
        params.(param_name) = value;
    end
    next_line = fgetl(fid);
end

fprintf('Loading options from %s... done.\n',opts_filename);
fclose(fid);
end