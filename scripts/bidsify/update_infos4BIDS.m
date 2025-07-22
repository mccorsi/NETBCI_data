%% Code to reorganize metadata
% Author: Marie-Constance Corsi
% Last modification: July, 22nd 2025

close all       % Close all open windows
clear all       % Clear all variables from the workspace
restoredefaultpath

%% Paths
root = pwd;
script_path = strcat(root,'/scripts/');
raw_path = strcat(root,'/NETBCI_db/'); % Folder with pre-BIDS data
bidsroot = strcat(root,'/NETBCI_postBIDS_db'); % Folder for output files
subj_data_path = strcat(raw_path, '/Metadata');
ftpath = strcat(root,'/fieldtrip-20240614'); % FieldTrip folder

addpath(ftpath)
ft_defaults 

%% Define subjects and sessions
% All subject relevant data is now stored in the 'metadata.mat' file. Import and use for the bidsify loop.
     
load(fullfile(subj_data_path, 'metadata'));
filenames = readtable(fullfile(bidsroot, 'filenames.csv'), 'Delimiter', ',' );

% Reorder according to anonymous id to be safe.
dates = filenames.acq_time;

%% General info 
n_sessions = 4;
n_subjects = 19;

% Common for all participants or info written to dataset_description file.
general = [];
general.method = 'copy';

% General metadata data
general.bidsroot                    = bidsroot;
general.datatype                    = 'meeg';
general.InstitutionName             = 'Paris Brain Institute';
general.InstitutionalDepartmentName = 'NERV Lab';
general.InstitutionAddress          = '47, boulevard de l hopital, 75013 Paris, France';

% Dataset description
general.dataset_description.Name            = 'NETBCI';
general.dataset_description.BIDSVersion     = 'v1.7.0';
general.dataset_description.EthicsApprovals = 'CPP-IDF-VI of Paris. 2016-A00626-45';
general.dataset_description.nsessions = n_sessions;
general.dataset_description.nsubjects = n_subjects;


% MEG recording info common to all recordings
general.meg.DewarPosition              = 'upright';    % REQUIRED. Position of the dewar during the MEG scan: "upright", "supine" or "degrees" of angle from vertical: for example on CTF systems, upright=15??, supine = 90??.
general.meg.SoftwareFilters = "SpatialCompensation_GradientOrder_0"; %, "GradientOrder", 0];       % REQUIRED. List of temporal and/or spatial software filters applied, orideally key:valuepairsofpre-appliedsoftwarefiltersandtheir parameter values: e.g., {"SSS": {"frame": "head", "badlimit": 7}}, {"SpatialCompensation": {"GradientOrder": Order of the gradient compensation}}. Write "n/a" if no software filters applied.
general.meg.DigitizedLandmarks         = true;         % REQUIRED. Boolean ("true" or "false") value indicating whether anatomical landmark points (i.e. fiducials) are contained within this recording.
general.meg.DigitizedHeadPoints        = true;         % REQUIRED. Boolean ("true" or "false") value indicating whether head points outlining the scalp/face surface are contained within this recording.
general.meg.RecordingType              = 'continuous';
general.meg.ContinuousHeadLocalization = false;
general.meg.PowerLineFrequency         = 50;           % REQUIRED. Frequency (in Hz) of the power grid at the geographical location of the MEG instrument (i.e. 50 or 60)
general.meg.Manufacturer = 'Elekta';
general.meg.SamplingFrequency = 250.0;
general.meg.DigitizedLandmarks = true;
general.meg.DigitizedHeadPoints = true;
general.meg.MEGChannelCount = 306;
general.meg.MEGREFChannelCount = 0;

% EEG recording info common to all recordings
general.eeg.SoftwareFilters = "n/a" ;     % REQUIRED. List of temporal and/or spatial software filters applied, orideally key:valuepairsofpre-appliedsoftwarefiltersandtheir parameter values: e.g., {"SSS": {"frame": "head", "badlimit": 7}}, {"SpatialCompensation": {"GradientOrder": Order of the gradient compensation}}. Write "n/a" if no software filters applied.
general.eeg.RecordingType              = 'continuous';
general.eeg.PowerLineFrequency         = 50;           % REQUIRED. Frequency (in Hz) of the power grid at the geographical location of the MEG instrument (i.e. 50 or 60)
general.eeg.Manufacturer = 'Elekta';
general.eeg.SamplingFrequency = 250.0;
general.eeg.EEGChannelCount = 74;

txt=jsonencode(general,PrettyPrint=true);
% remove white-spaces inside vectors and matrices
txt = regexprep(txt,',\s+(?=\d)',','); % , white-spaces digit
txt = regexprep(txt,',\s+(?=-)',','); % , white-spaces minussign
txt = regexprep(txt,'[\s+(?=\d)','['); % [ white-spaces digit
txt = regexprep(txt,'[\s+(?=-)','['); % [ white-spaces minussign
txt = regexprep(txt,'(?<=\d)\s+]',']'); % digit white-spaces ]

% Write the string to file
saveJSONfile(general, strcat(bidsroot,'/dataset_description.json'));

