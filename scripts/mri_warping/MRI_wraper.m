%% Warp template MRI to subject MRI for creating source model

% Adapted by M.-C. Corsi from Vinding, M. C., & Oostenveld, R. (2022). Sharing individualised template MRI data for MEG source reconstruction: A solution for open data while keeping subject confidentiality. NeuroImage, 119165. https://doi.org/10.1016/j.neuroimage.2022.119165
% Last modification: July 22nd, 2025

clear all; close all; clc; 

%% Paths
root = pwd;
script_path = strcat(root,'/scripts/');
raw_folder = strcat(root,'/NETBCI_db/'); % Folder with raw MRI data
out_folder = strcat(root,'/NETBCI_db/mri_wraped_tobidsify'); % Folder for output files
ftpath = strcat(root,'/fieldtrip-20240614'); % FieldTrip folder

addpath(ftpath)
ft_defaults 

%% Sessions - here, even though the MRI wad obtained at session 4 we have to align MRI with MEG coordinates for each session separately
sessions = {'Session1', 'Session2', 'Session3', 'Session4'};

%% Subjects
for kk_subj = 1:19 % subjects 1 & 2 done
   
   % Paths - specific to the subject (not the session since only one MRI has been acquired per subject)
    mri_path = fullfile(raw_folder, subj, 'MRI');  % Raw data folder
    sub_path = fullfile(out_folder);  % Output folder

   %% STEP 1A: Load subject MRI - do it once / subject
    % Load the subject anatomical image. Determine coordinate systen (ras, origin not a landmark).

    % Read MRI
    tmp_list = dir(fullfile(mri_path, '*.nii')); 
    filename_mri = tmp_list.name;
    raw_fpath = fullfile(mri_path, filename_mri); 
    mri_orig = ft_read_mri(raw_fpath);

    % Define coordinates of raw (r-a-s-n)
    mri_orig = ft_determine_coordsys(mri_orig, 'interactive', 'yes');

    %Save
    save(fullfile(sub_path, [subj, '_mri_orig.mat']), 'mri_orig')

   %% STEP 2: Convert subject MRI to desired coordinate system
    % % if needed, (Re)load data
    % load(fullfile(sub_path, [subj, '_mri_orig.mat']))

    %% Step 2A: Convert to Neuromag coordinate system
    cfg = [];
    cfg.method      = 'interactive';
    cfg.coordsys    = 'neuromag';
    mri_neuromag = ft_volumerealign(cfg, mri_orig);       

    % Save
    save(fullfile(sub_path, [subj, '_mri_neuromag.mat']), 'mri_neuromag')

end

%% specific to the session 
for kk_subj = 1:19 
    subj    = strcat('subject_',sprintf('%02d', kk_subj)); 

    mri_path = fullfile(raw_folder, subj, 'MRI');  % Raw data folder
    sub_path = fullfile(out_folder);  % Output folder

    load(fullfile(sub_path, [subj, '_mri_neuromag.mat']))
    load(fullfile(sub_path, [subj, '_mri_orig.mat']))


   %% Paths - specific to the subject (not the session)
    mri_path = fullfile(raw_folder, subj, 'MRI');  % Raw data folder
    sub_path = fullfile(out_folder);  % Output folder
    for kk_sess = 1:4  
        session = sessions{kk_sess};
        session_out = strcat('ses-',sprintf('%02d', kk_sess));
    
        %% Paths - specific to the session
        meg_path = fullfile(raw_folder, subj, session); % Raw data folder

        infname      = fullfile(out_folder, [subj, '_', session_out, 'mri_resliced.mat']);
        outfname_mat = fullfile(out_folder, [subj, '_', session_out, '_mri_warptmp.mat']);
        outfname_nii = fullfile(out_folder, [subj, '_', session_out, '_mri_warptmp']);
        
        
        %% Step 2B: Align MRI and MEG headpoints in MEG coordinate system (neuromag)
        % Get headshapes and sensor info from raw MEG file
        rawfile     = fullfile(meg_path, 'ds_testing06_trans_tsss.fif'); % Name of datafile with sensor and head points.
        headshape   = ft_read_headshape(rawfile);
        
        % Make sure units are mm
        headshape = ft_convert_units(headshape, 'mm');
        
        % Save headshapes and sensor info (optional)
        save(fullfile(meg_path, [subj, '-', session_out, '_headshape']), 'headshape')
        
        % Aligh MRI to MEG headpoints
        cfg = [];
        cfg.method              = 'headshape';
        cfg.headshape.headshape = headshape;
        cfg.headshape.icp       = 'yes';
        cfg.coordsys            = 'neuromag';
        mri_org_realign = ft_volumerealign(cfg, mri_neuromag);
        
        % Inspection (plot only)
        cfg.headshape.icp = 'no';
        ft_volumerealign(cfg, mri_org_realign);
        
        % Save
        save(fullfile(meg_path, [subj, '-', session_out, '_mri_org_realign']), 'mri_org_realign')
        disp('done')
        
        %% Step 2C: Reslice aligned image
        % Reslice to new coordinate system
        mri_org_resliced = ft_volumereslice([], mri_org_realign);
        
        % Save
        fprintf('saving...');
        save(fullfile(sub_path, [subj, '-', session_out, '_mri_org_resliced']), 'mri_org_resliced');
        disp('done')
        
        %% Step 3: Write subject volume as the "template". 

        % if necessary, (re)load data
        % load(fullfile(sub_path, [subj, '-', session_out, '_mri_org_resliced']));
    
        %% Step 3A: Write Neuromag image
        cfg = [];
        cfg.filetype    = 'nifti';          
        cfg.parameter   = 'anatomy';
        cfg.datatype    = 'double';
        cfg.filename    = fullfile(sub_path, [subj, '-', session_out, '_orig_neuromag_rs']);
        ft_volumewrite(cfg, mri_org_resliced)
    
        %% Step 3B: Make subject tissue probability maps
        cfg = [];
        cfg.spmmethod   = 'new';
        cfg.output      = 'tpm';
        cfg.write       = 'no';             
        org_seg = ft_volumesegment(cfg, mri_org_resliced); disp('done')
        
        % Rearrange data stucture for saving
        sub_tpm = mri_org_resliced;
        sub_tpm.anatomy          = org_seg.gray;
        sub_tpm.anatomy(:,:,:,2) = org_seg.white;
        sub_tpm.anatomy(:,:,:,3) = org_seg.csf;
        sub_tpm.anatomy(:,:,:,4) = org_seg.bone;
        sub_tpm.anatomy(:,:,:,5) = org_seg.softtissue;
        sub_tpm.anatomy(:,:,:,6) = org_seg.air;
        sub_tpm.dim = [org_seg.dim, size(sub_tpm.anatomy, 4)];
    
        %% Step 3C: Write TPM volume
        cfg = [];
        cfg.filetype    = 'nifti';          
        cfg.parameter   = 'anatomy';
        cfg.datatype    = 'double';
        cfg.filename    = fullfile(sub_path, [subj, '-', session_out, '_sub_tpm']);
        ft_volumewrite(cfg, sub_tpm);
    
        %% STEP 4: warp a template MRI to the individual MRI "template"
        % Load template MRI
        load standard_mri           % Load Colin 27 from FieldTrip
        mri_colin = mri;            % Rename to avoid confusion
    
        %% Step 3A: Do initial alignmet of fiducials to target coordsys (optional)
        cfg = [];
        cfg.method      = 'interactive';
        cfg.coordsys    = 'neuromag';
        mri_colin_neuromag = ft_volumerealign(cfg, mri_colin);     
    
        %% Step 3B: Normalise template -> subject (neuromag coordsys)
        cfg = [];
        cfg.nonlinear        = 'yes';       % Non-linear warping
        cfg.spmmethod        = 'new';       % Note: method="new" uses SPM's default posterior tissue maps unless we specify 'cfg.tpm'.
        cfg.spmversion       = 'spm12';     % Default = "spm12"
        cfg.templatecoordsys = 'neuromag';  % Coordinate system of the template
        cfg.template         = fullfile(sub_path, [subj, '-', session_out, '_orig_neuromag_rs.nii']); % Subject "template" created in step 3A
        cfg.tpm              = fullfile(sub_path, [subj, '-', session_out, '_sub_tpm.nii']);         % Subject TPM created in step 3C
        mri_warptmp = ft_volumenormalise(cfg, mri_colin_neuromag);
        
        % Determine unit of volume (mm)
        mri_warptmp = ft_determine_units(mri_warptmp);
        
        % Plot for inspection
        ft_sourceplot([], mri_warptmp); title('Warped2neuromag')
    
        %% Save
        fprintf('saving...')
        save(fullfile(sub_path, [subj, '-', session_out, '_mri_warptmp']), 'mri_warptmp')
        disp('done')
    
        %% Write image
        cfg = [];
        cfg.filetype    = 'nifti';          
        cfg.parameter   = 'anatomy';
        cfg.datatype    = 'double';
        cfg.filename    = outfname_nii;
        ft_volumewrite(cfg, mri_warptmp)
        disp('done')
        
        fprintf('%s done\n', subj)
        clear mri_warptmp sub_tmp org_seg mri_resliced
        close all
    end
    clear mri_orig mri_neuromag
    disp(strcat(subj, ': done (all sessions)'))
end
