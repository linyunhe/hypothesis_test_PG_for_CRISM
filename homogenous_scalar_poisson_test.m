% Codes for CRISM L data only now
% Linyun He
% 03/13/2019
function res = homogenous_scalar_poisson_test()
    cur_folder = pwd;
    [ssa_filename,ssa_path] = uigetfile({'*.img*','*.bsq*'},'SSA file');
    [wa_filename,wa_path] = uigetfile({'*.img*','*.bsq*'},'WA file');
    output_folder = uigetdir(cur_folder,'Output folder');
    video_filename = fullfile(output_folder,'Poisson_test_plot_cdfpdf.avi');
    measure = multibandread_fromENVIHeader(fullfile(ssa_path,ssa_filename));
    wa = read_envi_data(fullfile(wa_path,wa_filename));
    
    band_min = 188;
    band_max = 433;
    [~,num_col,num_band] = size(measure);
    if num_col == 640
        col_min = 32;
        col_max = 631;
    else
        col_min = 1;
        col_max = num_col;
    end
    measure = measure(:,col_min:col_max,:);
    measure(measure>1)=0;
    wa_use = squeeze(wa(1,round(num_band/2),band_min:band_max))/1000;
    [mask_new,scalar,spat_fig,scalar_fig] = cutpart(measure,[24,174,242],1,wa_use);
    
    samples_num = sum(mask_new(:))/num_band;
    region = reshape(measure(mask_new),[samples_num,num_band]);
    scaled_region = region.*scalar;
    savefig(scalar_fig,fullfile(output_folder,'Scalars_bands.fig'));
    savefig(spat_fig,fullfile(output_folder,'Spat_area_homo.fig'));
    
    writerobj = VideoWriter(video_filename);
    writerobj.FrameRate = 2;
    open(writerobj);
    hfig = figure;
    for band = 1:num_band
        tmp = scaled_region(:,band);
        tmp = tmp(tmp > 0);
        [nn,xx] = hist(tmp,30);
        nncdf = cumsum(nn)./sum(nn);
        nnpdf = nn./sum(nn);
        space = mean(diff(xx));
        poisshat = mean(tmp);
        ypoisscdf = poisscdf(xx+space/2,poisshat);
        ypoisspdf = normpdf(xx,poisshat,sqrt(poisshat))*space;
        tmpcdf = zeros(2,30);
        tmppdf = tmpcdf;
        tmpcdf(1,:) = nncdf;
        tmpcdf(2,:) = ypoisscdf;
        tmppdf(1,:) = nnpdf;
        tmppdf(2,:) = ypoisspdf;
        % figure;
        [ax,h1,h2] = plotyy(xx,tmpcdf,xx,tmppdf);
        set(ax(1),'ycolor',[0 0 0.5]);
        set(ax(2),'ycolor',[0 0.5 0]);
        set(h1(1),'color',[0 0 0.5],'linestyle','none','marker','d');
        set(h1(2),'color',[0 0 0.5],'linestyle','-','marker','none');
        set(h2(2),'color',[0 0.5 0],'linestyle','-','marker','none');
        set(h2(1),'color',[0 0.5 0],'linestyle','none','marker','x');
        set(get(ax(1),'YLabel'),'Str','Cumulative Distribution Function');
        set(get(ax(2),'YLabel'),'Str','Probability Density Function');
        set(ax(1),'ytick',0:.1:1);
        set(ax(2),'ytick',0:.02:0.1);
        set(ax(2),'YLim',[0,0.11]);
        xlabel('SSA/\alpha');
        title(sprintf('Distribution of homogenous area when band is %i',band));
        legend('Measurement CDF','Theoritical CDF','Measurement PDF','Theoritical PDF','Location','Best');
        frame = getframe(hfig);
        writeVideo(writerobj,frame);
    end
    close(writerobj);
    res = true;
end