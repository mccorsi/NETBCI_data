%% Code to reformat files that contain events
% Author: Marie-Constance Corsi
% Last modification: July, 22nd 2025

root = pwd;
script_path = strcat(root,'/scripts/');
path_data2bidsify_root = strcat(root,'/NETBCI_db/'); 
path_dataPostbidsify_root = strcat(root,'/NETBCI_postBIDS_db'); 

% MI: motor imagery trials, BASELINE: resting state trials
list_events_id = {'TRIAL_MI', 'TRIAL_BASELINE'};

%% Batch across subjects and sessions
for kk_subj = 1:19
    subj_id = sprintf( '%02d', kk_subj );
    for kk_session = 1:4
        
        path2save = strcat(path_data2bidsify_root, 'subject_' , subj_id, '/Session', num2str(kk_session), '/');

        for kk_run = 1:6
            run_id = sprintf( '%02d', kk_run );
            try
                filename2load = strcat(path2save, 'ft_events_ds_testing', run_id, '_trans_tsss_checked');
                DoMyReformat_events_ft2MNE(filename2load, list_events_id)
            catch
                filename2load = strcat(path2save, 'ft_events_ds_testing', run_id, '_tsss_checked');
                DoMyReformat_events_ft2MNE(filename2load, list_events_id)
            end
            
        end
    end
end



