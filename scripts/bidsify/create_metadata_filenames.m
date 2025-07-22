%% Code to create metadata with relevant information from a personal .mat file
% Author: Marie-Constance Corsi
% Last modification: July, 22nd 2025

clear all; clc;


root = pwd;
script_path = strcat(root,'/scripts/');
raw_path = strcat(root,'/NETBCI_db/'); % Folder with pre-BIDS data
bidsroot = strcat(root,'/NETBCI_postBIDS_db'); % Folder for output files
metadata_path = strcat(raw_path, '/Metadata');
ftpath = strcat(root,'/fieldtrip-20240614'); % FieldTrip folder


load('/Users/marieconstance.corsi/Documents/Personal/behavior.mat'); % personal file that contains information
load(strcat(metadata_path,'GeneralInfos.mat')); % personal file that contains information

metadata.sex = sex([1:15,17:20]); 
metadata.agebin = age ([1:15,17:20]);

metadata.BCIScoresPerRun.session1 = behavior_updated.BCI.Perf.Sess1.Runs([1:15,17:20]);
metadata.BCIScoresPerRun.session2 = behavior_updated.BCI.Perf.Sess2.Runs([1:15,17:20]);
metadata.BCIScoresPerRun.session3 = behavior_updated.BCI.Perf.Sess3.Runs([1:15,17:20]);
metadata.BCIScoresPerRun.session4 = behavior_updated.BCI.Perf.Sess4.Runs([1:15,17:20]);
metadata.rosenbergScale = behavior_updated.PsyQuest.RosenbergScale([1:15,17:20]);


metadata.EMG28.amotivation = behavior_updated.PsyQuest.EMG28.Amotivation([1:15,17:20]);
metadata.EMG28.extrinsic.externReg = behavior_updated.PsyQuest.EMG28.Extrinsic.externReg([1:15,17:20]);
metadata.EMG28.extrinsic.identifReg = behavior_updated.PsyQuest.EMG28.Extrinsic.identifReg([1:15,17:20]);
metadata.EMG28.extrinsic.introjReg = behavior_updated.PsyQuest.EMG28.Extrinsic.introjReg([1:15,17:20]);
metadata.EMG28.intrinsic.accomplishment = behavior_updated.PsyQuest.EMG28.Intrinsic.accomplishment([1:15,17:20]);
metadata.EMG28.intrinsic.knowledge = behavior_updated.PsyQuest.EMG28.Intrinsic.knowledge([1:15,17:20]);
metadata.EMG28.intrinsic.stimulation = behavior_updated.PsyQuest.EMG28.Intrinsic.stimulation([1:15,17:20]);


metadata.VMIQ2.externalVisualImagery = behavior_updated.PsyQuest.VMIQ_2.ExternalVisualImagery([1:15,17:20]);
metadata.VMIQ2.InternalVisualImagery = behavior_updated.PsyQuest.VMIQ_2.InternalVisualImagery([1:15,17:20]);
metadata.VMIQ2.KinaestheticImagery = behavior_updated.PsyQuest.VMIQ_2.KinaestheticImagery([1:15,17:20]);
metadata.STAI_YA.session1 = behavior_updated.PsyQuest.STAI_YA.Sess1([1:15,17:20]);
metadata.STAI_YA.session2 = behavior_updated.PsyQuest.STAI_YA.Sess2([1:15,17:20]);
metadata.STAI_YA.session3 = behavior_updated.PsyQuest.STAI_YA.Sess3([1:15,17:20]);
metadata.STAI_YA.session4 = behavior_updated.PsyQuest.STAI_YA.Sess4([1:15,17:20]);

json_data = jsonencode(metadata);
writematrix(json_data, strcat(metadata_path,"metadata.csv"));
writematrix(json_data, strcat(metadata_path,"metadata.txt"));

save(strcat(metadata_path,"metadata.mat"), "metadata")

writetable(struct2table(metadata), 'Structure_Example.csv')

%% create filenames.csv from the scans -> cf 'Update_Participants_Experiments_Infos.py'
