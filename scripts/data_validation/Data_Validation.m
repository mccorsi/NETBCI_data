%% Code to validate the data
% Author: Marie-Constance Corsi
% Last modification: July, 22nd 2025
clear all
addpath('/Users/marieconstance.corsi/Documents/GitHub/NETBCI_data/packages/fieldtrip-20240614/')
ft_defaults

%% Paths
root = pwd;
script_path = strcat(root,'/scripts/');
bids_path = strcat(root,'/NETBCI_db/');
fig_path = strcat(root,'/figures/');

addpath(genpath(script_path)); 
%% Load subjects (not mandatory)
T = readtable(fullfile(bids_path, 'linkdata.csv'), 'VariableNamingRule', 'preserve');
subjects = T.anonymized_id;

%% Compute PSD 
for kk_sess = 1:4
    sess = strcat('ses-0', num2str(kk_sess));
    for kk_subj = 1:length(subjects)
        subj = [subjects{kk_subj}];
        disp(strcat('Processing :  ', subjects{kk_subj}, '-', sess));
        subjdir = fullfile(bids_path, subj, sess, 'meg');

        infile = dir(strcat(subjdir, '/*task-rest*','.fif'));
        list_filenames = {infile.name};

        for kk_file = 1:length(infile)
            cfg = [];
            cfg.dataset = fullfile(subjdir, list_filenames{kk_file});
            cfg.channel = 'MEG';
            cfg.checkmaxfilter = 'no'; 
            raw = ft_preprocessing(cfg);
        
            cfg = [];
            cfg.length = 2;
            cfg.overlap = 0.5;
            epo = ft_redefinetrial(cfg, raw);
        
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.pad = 'nextpow2';
            cfg.output = 'pow';
            cfg.taper = 'hanning';
            cfg.foilim = [0.5, 100];
            freq = ft_freqanalysis(cfg, epo);
        
            allFreq{kk_sess, kk_subj, kk_file} = freq;
        
            cfg = [];
            cfg.avgoverfreq = 'no';
            cfg.avgoverchan = 'yes';
            freqGlob = ft_selectdata(cfg, freq);
        
            allFreqGlob(kk_sess, kk_subj, kk_file,:) = freqGlob.powspctrm;
        end
    end
    
end
disp('DONE')
freqax = freqGlob.freq;

%% Save
fprintf('Saving... ')
save(strcat(root,'/misc/freq.mat'), 'allFreq', 'allFreqGlob', 'freqax', '-v7.3')
disp('done')
 
%% Average plot per session and per resting state file (ie beginning and end of the session)
for kk_sess = 1:4
    avgPsd_RS1(kk_sess,:) = squeeze(mean(allFreqGlob(kk_sess, :, 1,:),2));
    avgPsd_RS2(kk_sess,:) = squeeze(mean(allFreqGlob(kk_sess, :, 2,:),2));
end

fig = figure(1);
for kk_sess = 1:4
     
    if kk_sess == 1
        col = '#665656';
    elseif kk_sess == 2 
        col = '#BA756C';
    elseif kk_sess == 3
        col = '#C7AA99';
    else
        col = '#DAC188';
    end
    
    subplot(2,4,kk_sess)
    gr1 = plot(freqax, log(squeeze(avgPsd_RS1(kk_sess,:))), 'color', col, 'LineWidth', 3.6); hold on
    for kk_subj = 1:size(allFreqGlob,2)
        plot(freqax, log(squeeze(allFreqGlob(kk_sess, kk_subj,1,:))), 'color', col) % RS1
        hold on;
    end
    xlim([0,60]); ylim([-60,-52.5])
    xlabel('Frequency (Hz)'); ylabel('Log-Power')
    title(strcat('Session #', num2str(kk_sess)))
    ax = gca;
    ax.FontWeight = 'bold';
    box off


    subplot(2,4,kk_sess+4)
    gr1 = plot(freqax, log(squeeze(avgPsd_RS2(kk_sess,:))), 'color', col, 'LineWidth', 3.6); hold on
    for kk_subj = 1:size(allFreqGlob,2)
        plot(freqax, log(squeeze(allFreqGlob(kk_sess, kk_subj,2,:))), 'color', col) % RS1
        hold on;
    end
    % title(strcat('Post-experiment - Session #', num2str(kk_sess)))
    xlim([0,60]); ylim([-60,-52.5])
    xlabel('Frequency (Hz)'); ylabel('Log-Power')
    ax = gca;
    ax.FontWeight = 'bold';
    box off
end
% Add row titles
annotation('textbox', [0.01, 0.68, 0.1, 0.1], 'String', 'Pre-experiment', 'EdgeColor', 'none', 'FontSize', 13);%, 'FontWeight', 'bold');
annotation('textbox', [0.01, 0.20, 0.1, 0.1], 'String', 'Post-experiment', 'EdgeColor', 'none', 'FontSize', 13);%, 'FontWeight', 'bold');


% Save
savefig(fig, fullfile(fig_path, 'psd.fig'));
saveas(fig, fullfile(fig_path, 'psd.tif'), 'tiff')
exportgraphics(fig, fullfile(fig_path, 'psd.tif'), 'Resolution',600)
disp('done)')

%% Topoplots
for kk_sess=1:4
    % Pre-experiment
    cfg = [];
    avgTopo_RS1{kk_sess} = ft_freqgrandaverage(cfg, allFreq{kk_sess, :, 1});
    
    fidx = find(avgTopo_RS1{kk_sess}.freq==10);
    ul = max([avgTopo_RS1{kk_sess}.powspctrm(:, fidx); avgTopo_RS1{kk_sess}.powspctrm(:, fidx)]);
    ll = min([avgTopo_RS1{kk_sess}.powspctrm(:, fidx); avgTopo_RS1{kk_sess}.powspctrm(:, fidx)]);
    zlim10 = [ll, ul];
    
    fidx = find(avgTopo_RS1{kk_sess}.freq==20);
    ul = max([avgTopo_RS1{kk_sess}.powspctrm(:, fidx); avgTopo_RS1{kk_sess}.powspctrm(:, fidx)]);
    ll = min([avgTopo_RS1{kk_sess}.powspctrm(:, fidx); avgTopo_RS1{kk_sess}.powspctrm(:, fidx)]);
    zlim20 = [ll, ul];

    % Post-experiment
    cfg = [];
    avgTopo_RS2{kk_sess} = ft_freqgrandaverage(cfg, allFreq{kk_sess, :, 2});
    
    fidx = find(avgTopo_RS2{kk_sess}.freq==10);
    ul = max([avgTopo_RS2{kk_sess}.powspctrm(:, fidx); avgTopo_RS2{kk_sess}.powspctrm(:, fidx)]);
    ll = min([avgTopo_RS2{kk_sess}.powspctrm(:, fidx); avgTopo_RS2{kk_sess}.powspctrm(:, fidx)]);
    zlim10 = [ll, ul];
    
    fidx = find(avgTopo_RS2{kk_sess}.freq==20);
    ul = max([avgTopo_RS2{kk_sess}.powspctrm(:, fidx); avgTopo_RS2{kk_sess}.powspctrm(:, fidx)]);
    ll = min([avgTopo_RS2{kk_sess}.powspctrm(:, fidx); avgTopo_RS2{kk_sess}.powspctrm(:, fidx)]);
    zlim20 = [ll, ul];

end


cfg = [];
cfg.parameter   = 'powspctrm';
cfg.xlim        = [8, 13];
cfg.layout      = 'neuromag306mag.lay';
cfg.comment     = 'no';
cfg.colormap =   ft_colormap('viridis'); %('-RdBu');
cfg.zlim        = [0, 7e-27];
cfg.figure      = 'no';


for kk_sess=1:4
    
    fig_pre = figure();
    ft_topoplotER(cfg, avgTopo_RS1{kk_sess}); colorbar
    exportgraphics(fig_pre,strcat(fig_path,'topo_alpha_pre_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

    fig_post = figure();
    ft_topoplotER(cfg, avgTopo_RS2{kk_sess}); colorbar
    exportgraphics(fig_post,strcat(fig_path,'topo_alpha_post_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

end

%% BONUS - topoplots in other frequency bands
cfg = [];
cfg.parameter   = 'powspctrm';
cfg.xlim        = [8, 10];
cfg.layout      = 'neuromag306mag.lay';
cfg.comment     = 'no';
cfg.colormap =   ft_colormap('viridis'); %('-RdBu');
cfg.zlim        = [0, 9e-27];
cfg.figure      = 'no';


for kk_sess=1:4
    
    fig_pre = figure();
    ft_topoplotER(cfg, avgTopo_RS1{kk_sess}); colorbar
    exportgraphics(fig_pre,strcat(fig_path,'topo10_pre_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

    fig_post = figure();
    ft_topoplotER(cfg, avgTopo_RS2{kk_sess}); colorbar
    exportgraphics(fig_post,strcat(fig_path,'topo10_post_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

end

cfg = [];
cfg.parameter   = 'powspctrm';
cfg.xlim        = [15, 20];
cfg.layout      = 'neuromag306mag.lay';
cfg.comment     = 'no';
cfg.colormap =   ft_colormap('viridis');
cfg.zlim        = [0, 9e-28];
cfg.figure      = 'no';

for kk_sess=1:4
    
    fig_pre = figure();
    ft_topoplotER(cfg, avgTopo_RS1{kk_sess}); colorbar
    exportgraphics(fig_pre,strcat(fig_path,'topo20_pre_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

    fig_post = figure();
    ft_topoplotER(cfg, avgTopo_RS2{kk_sess}); colorbar
    exportgraphics(fig_post,strcat(fig_path,'topo20_post_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

end


cfg = [];
cfg.parameter   = 'powspctrm';
cfg.xlim        = [30, 40];
cfg.layout      = 'neuromag306mag.lay';
cfg.comment     = 'no';
cfg.colormap =   ft_colormap('viridis');
cfg.zlim        = [5e-29, 12e-29];
cfg.figure      = 'no';

for kk_sess=1:4
    
    fig_pre = figure();
    ft_topoplotER(cfg, avgTopo_RS1{kk_sess}); colorbar
    exportgraphics(fig_pre,strcat(fig_path,'topo40_pre_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

    fig_post = figure();
    ft_topoplotER(cfg, avgTopo_RS2{kk_sess}); colorbar
    exportgraphics(fig_post,strcat(fig_path,'topo40_post_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

end


cfg = [];
cfg.parameter   = 'powspctrm';
cfg.xlim        = [50, 50];
cfg.layout      = 'neuromag306mag.lay';
cfg.comment     = 'no';
cfg.colormap =   ft_colormap('viridis');
cfg.zlim        = [5e-29, 21e-29];
cfg.figure      = 'no';

for kk_sess=1:4
    
    fig_pre = figure();
    ft_topoplotER(cfg, avgTopo_RS1{kk_sess}); colorbar
    exportgraphics(fig_pre,strcat(fig_path,'topo50_pre_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

    fig_post = figure();
    ft_topoplotER(cfg, avgTopo_RS2{kk_sess}); colorbar
    exportgraphics(fig_post,strcat(fig_path,'topo50_post_sess-',num2str(kk_sess),'.png'), 'Resolution',600)

end
% % END
