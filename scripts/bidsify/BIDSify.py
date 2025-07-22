##% Code to bidsify M-EEG data
# Author: Marie-Constance Corsi
# Last modification: July, 22nd 2025

import os
import os.path as op
import shutil

# %%
import mne
import numpy as np
import pandas as pd
import scipy.io as sio
from mne_bids import (
    BIDSPath,
    make_report,
    write_raw_bids,
)
from mne_bids.stats import count_events


# %% Functions
def _convert_events_ft2mne(
    sample_data_raw_file, filename_events_mat_MI, filename_events_mat_Rest
):
    # load first file of events for the first type
    events_fname_MI = op.join(filename_events_mat_MI)
    tmp_events_contents_MI = sio.loadmat(events_fname_MI)
    t_events_contents_MI = np.array(tmp_events_contents_MI["mat_events"])
    col_id_MI = np.reshape(
        [1] * len(t_events_contents_MI), (len(t_events_contents_MI), 1)
    )
    events_contents_MI = np.hstack((t_events_contents_MI, col_id_MI))

    events_fname_Rest = op.join(filename_events_mat_Rest)
    tmp_events_contents_Rest = sio.loadmat(events_fname_Rest)
    t_events_contents_Rest = np.array(tmp_events_contents_Rest["mat_events"])
    col_id_Rest = np.reshape(
        [2] * len(t_events_contents_Rest), (len(t_events_contents_Rest), 1)
    )
    events_contents_Rest = np.hstack((t_events_contents_Rest, col_id_Rest))

    tmp_events = np.concatenate((events_contents_MI, events_contents_Rest))
    events_array = tmp_events[tmp_events[:, 0].argsort()]

    filename_event_file = sample_data_raw_file + "-eve.fif"
    mne.write_events(filename_event_file, events_array, overwrite=True, verbose=None)
    print("Saving the MNE event file...")

    return filename_event_file, events_array


# %% Paths
if os.path.basename(os.getcwd()) == "scripts":
    os.chdir("..")
basedir = os.getcwd()

os.chdir(basedir)
path_data2bidsify_root = basedir + "/NETBCI_db/"
path_dataPostbidsify_root = basedir + "/NETBCI_postBIDS_db/"
list_mod = ["MEG", "EEG"]

# %% General infos
ch_names_eeg = [
    "FP1",
    "FPz",
    "FP2",
    "AF1",
    "AF3",
    "AFz",
    "AF4",
    "AF8",
    "F7",
    "F5",
    "F3",
    "F1",
    "Fz",
    "F2",
    "F4",
    "F6",
    "F8",
    "FT9",
    "FT7",
    "FC5",
    "FC3",
    "FC1",
    "FCz",
    "FC2",
    "FC4",
    "FC6",
    "FT8",
    "FT10",
    "T9",
    "T7",
    "C5",
    "C3",
    "C1",
    "Cz",
    "C2",
    "C4",
    "C6",
    "T8",
    "T10",
    "TP9",
    "TP7",
    "CP5",
    "CP3",
    "CP1",
    "CPz",
    "CP2",
    "CP4",
    "CP6",
    "TP8",
    "TP10",
    "P9",
    "P7",
    "P3",
    "P1",
    "Pz",
    "P2",
    "P4",
    "P6",
    "P8",
    "P10",
    "PO9",
    "PO7",
    "PO5",
    "PO3",
    "POz",
    "PO4",
    "PO8",
    "PO10",
    "O1",
    "Oz",
    "O2",
    "O9",
    "IZ",
    "O10",
]

participantsInfos = sio.loadmat(path_data2bidsify_root + "GeneralInfos.mat")
age = participantsInfos["age"][0]
sex = participantsInfos["sex"][0]

power_line = 50

event_dict = {"MI": 1, "Rest": 2}
# metadata for each epoch shall include events from the range: [0.0, 5] s (presentation of the target + feedback, t=-1:0 ISI, and t : [5, 7] s results of the trial
metadata_tmin, metadata_tmax = 0.0, 5.0

bids_root = path_dataPostbidsify_root

list_subj = list(range(1, 20))

# %% Convertion into bids
if op.exists(bids_root):
    shutil.rmtree(bids_root)

raw_eeg_list = list()
bids_eeg_list = list()
raw_meg_list = list()
bids_meg_list = list()

files2debug = pd.DataFrame(
    columns=["Subject", "Session", "Modality", "Recording", "Filename"]
)

for kk_subj in list_subj:
    subject_id = "{:02d}".format(
        kk_subj
    )  # zero padding to account for >10 subjects in this dataset
    tmp_sex = int(sex[kk_subj - 1])
    tmp_age = int(age[kk_subj - 1])

    for kk_session in list(range(1, 5)):
        session_id = "{:01d}".format(kk_session)
        path_subj_sess = (
            path_data2bidsify_root + "subject_" + subject_id + "/Session" + session_id
        )

        task = "rest"
        for kk_rest in [1, 2]:
            rs_id = "{:02d}".format(kk_rest)
            sample_data_raw_file = (
                path_subj_sess + "/ds_resting_state" + rs_id + "_trans_tsss"
            )
            try:
                raw = mne.io.read_raw_fif(sample_data_raw_file + ".fif")
                raw.info["line_freq"] = (
                    power_line  # specify power line frequency as required by BIDS
                )
                subject_info = raw.info.get("subject_info", None)
                subject_info["age"] = tmp_age
                subject_info["sex"] = tmp_sex
                raw.info["subject_info"] = subject_info
                # check dowsampling
                if raw.info["sfreq"] != 250:
                    raw.resample(sfreq=250)
                raw_eeg = raw.copy()
                raw_meg = raw.copy()

                # Get the electrode coordinates
                raw_eeg.rename_channels(
                    {
                        f"EEG{ii:03d}": ch_names_eeg_tmp
                        for ii, ch_names_eeg_tmp in enumerate(ch_names_eeg, 1)
                    }
                )
                mon = mne.channels.read_dig_fif(sample_data_raw_file + ".fif")
                mon.ch_names[1 : len(mon.ch_names)] = ch_names_eeg
                trans = mne.channels.compute_native_head_t(mon)
                raw_eeg.set_montage(mon)
                raw_eeg.pick(picks=["eeg"], exclude=["emg", "eog", "meg", "ecg"])

                raw_meg.pick(picks=["meg"], exclude=["emg", "eog", "eeg", "ecg"])
                # save results
                if kk_subj > 16:  # to avoid pb with subj 16 without shared data
                    subject_id_bids = "{:02d}".format(
                        kk_subj - 1
                    )  # zero padding to account for >10 subjects in this dataset
                else:
                    subject_id_bids = "{:02d}".format(kk_subj)
                session_id = "{:02d}".format(kk_session)
                raw_eeg_list.append(raw)
                bids_path_eeg = BIDSPath(
                    subject=subject_id_bids,
                    session=session_id,
                    task=task,
                    run=rs_id,
                    root=bids_root,
                    datatype="eeg",
                )
                bids_eeg_list.append(bids_path_eeg)

                raw_meg_list.append(raw_meg)
                bids_path_meg = BIDSPath(
                    subject=subject_id_bids,
                    session=session_id,
                    task=task,
                    run=rs_id,
                    root=bids_root,
                    datatype="meg",
                )
                bids_meg_list.append(bids_path_meg)
                try:
                    write_raw_bids(
                        raw=raw_eeg,
                        bids_path=bids_path_eeg,
                        montage=mon,
                        overwrite=True,
                        format="auto",
                    )
                    write_raw_bids(
                        raw=raw_meg,
                        bids_path=bids_path_meg,
                        overwrite=True,
                        format="auto",
                    )
                except:
                    print(
                        "Failed to write raw bids - subj:"
                        + str(kk_subj)
                        + "- sess: "
                        + str(kk_session)
                    )
                    data = {
                        "Subject": kk_subj,
                        "Session": kk_session,
                        "Modality": "MEEG",
                        "Recording": "RS" + str(rs_id),
                        "Filename": sample_data_raw_file,
                    }
                    data = pd.DataFrame(data.items())
                    data = data.transpose()
                    data.columns = data.iloc[0]
                    data = data.drop(data.index[[0]])
                    files2debug = pd.concat((files2debug, data))

            except:
                for kk_rest in [1, 2]:
                    rs_id = "{:02d}".format(kk_rest)
                    for kk_mod in [0, 1]:
                        modality_id = list_mod[kk_mod]
                        sample_data_raw_file = (
                            path_subj_sess
                            + "/ds_resting_state"
                            + rs_id
                            + "_trans_tsss_"
                            + modality_id
                            + ".fif"
                        )

                        if kk_mod == "EEG":
                            raw_eeg = mne.io.read_raw_fif(sample_data_raw_file + ".fif")
                            raw_eeg.info["line_freq"] = (
                                power_line  # specify power line frequency as required by BIDS
                            )
                            subject_info = raw_eeg.info.get("subject_info", None)
                            subject_info["age"] = tmp_age
                            subject_info["sex"] = tmp_sex
                            raw_eeg.info["subject_info"] = subject_info
                            # check dowsampling
                            if raw_eeg.info["sfreq"] != 250:
                                raw_eeg.resample(sfreq=250)
                        elif kk_mod == "MEG":
                            raw_meg = mne.io.read_raw_fif(sample_data_raw_file + ".fif")
                            raw_meg.info["line_freq"] = (
                                power_line  # specify power line frequency as required by BIDS
                            )
                            subject_info = raw_meg.info.get("subject_info", None)
                            subject_info["age"] = tmp_age
                            subject_info["sex"] = tmp_sex
                            raw_meg.info["subject_info"] = subject_info
                            # check dowsampling
                            if raw_meg.info["sfreq"] != 250:
                                raw_meg.resample(sfreq=250)

                    # save results
                    if kk_subj > 16:  # to avoid pb with subj 16 without shared data
                        subject_id_bids = "{:02d}".format(
                            kk_subj - 1
                        )  # zero padding to account for >10 subjects in this dataset
                    else:
                        subject_id_bids = "{:02d}".format(kk_subj)

                    session_id = "{:02d}".format(kk_session)
                    raw_eeg_list.append(raw)
                    bids_path_eeg = BIDSPath(
                        subject=subject_id_bids,
                        session=session_id,
                        task=task,
                        run=rs_id,
                        root=bids_root,
                        datatype="eeg",
                    )
                    bids_eeg_list.append(bids_path_eeg)

                    raw_meg_list.append(raw_meg)
                    bids_path_meg = BIDSPath(
                        subject=subject_id_bids,
                        session=session_id,
                        task=task,
                        run=rs_id,
                        root=bids_root,
                        datatype="meg",
                    )
                    bids_meg_list.append(bids_path_meg)

                    try:
                        write_raw_bids(
                            raw=raw_eeg,
                            bids_path=bids_path_eeg,
                            montage=mon,
                            overwrite=True,
                            format="auto",
                        )  # fif format does not work here...
                        write_raw_bids(
                            raw=raw_meg,
                            bids_path=bids_path_meg,
                            overwrite=True,
                            format="auto",
                        )
                    except:
                        print(
                            "Failed to write raw bids - subj:"
                            + str(kk_subj)
                            + "- sess: "
                            + str(kk_session)
                        )
                        data = {
                            "Subject": kk_subj,
                            "Session": kk_session,
                            "Modality": "MEEG",
                            "Recording": "RS" + str(rs_id),
                            "Filename": sample_data_raw_file,
                        }
                        data = pd.DataFrame(data.items())
                        data = data.transpose()
                        data.columns = data.iloc[0]
                        data = data.drop(data.index[[0]])
                        files2debug = pd.concat((files2debug, data))

        task = "MotorImageryRest"
        for kk_run in list(range(1, 7)):
            run_id = "{:02d}".format(kk_run)

            sample_data_raw_file = (
                path_subj_sess + "/ds_testing" + run_id + "_trans_tsss"
            )
            raw = mne.io.read_raw_fif(sample_data_raw_file + ".fif")
            raw.info["line_freq"] = (
                power_line  # specify power line frequency as required by BIDS
            )
            subject_info = raw.info.get("subject_info", None)
            subject_info["age"] = tmp_age
            subject_info["sex"] = tmp_sex
            raw.info["subject_info"] = subject_info

            raw_eeg = raw.copy()
            raw_meg = raw.copy()

            # load first file of events for the first type
            if os.path.isfile(sample_data_raw_file + "-eve.fif") == False:
                try:
                    events_fname_MI = op.join(
                        path_subj_sess
                        + "/ft_events_ds_testing"
                        + run_id
                        + "_trans_tsss_checked_TRIAL_MI.mat"
                    )
                    events_fname_Rest = op.join(
                        path_subj_sess
                        + "/ft_events_ds_testing"
                        + run_id
                        + "_trans_tsss_checked_TRIAL_BASELINE.mat"
                    )
                    filename_event_file, events_array = _convert_events_ft2mne(
                        sample_data_raw_file, events_fname_MI, events_fname_Rest
                    )
                except:
                    events_fname_MI = op.join(
                        path_subj_sess
                        + "/ft_events_ds_testing"
                        + run_id
                        + "_tsss_checked_TRIAL_MI.mat"
                    )
                    events_fname_Rest = op.join(
                        path_subj_sess
                        + "/ft_events_ds_testing"
                        + run_id
                        + "_tsss_checked_TRIAL_BASELINE.mat"
                    )
                    filename_event_file, events_array = _convert_events_ft2mne(
                        sample_data_raw_file, events_fname_MI, events_fname_Rest
                    )
            else:
                events_array = mne.read_events(
                    filename=sample_data_raw_file + "-eve.fif"
                )

            # auto-create metadata:
            metadata, events, event_id = mne.epochs.make_metadata(
                events=events_array,
                event_id=event_dict,
                tmin=metadata_tmin,
                tmax=metadata_tmax,
                sfreq=raw.info["sfreq"],
            )
            metadata_description = {c: i for i, c in enumerate(metadata.columns)}

            # Get the electrode coordinates
            raw_eeg.rename_channels(
                {
                    f"EEG{ii:03d}": ch_names_eeg_tmp
                    for ii, ch_names_eeg_tmp in enumerate(ch_names_eeg, 1)
                }
            )
            mon = mne.channels.read_dig_fif(sample_data_raw_file + ".fif")
            mon.ch_names[1 : len(mon.ch_names)] = ch_names_eeg
            trans = mne.channels.compute_native_head_t(mon)
            raw_eeg.set_montage(mon)

            # save results
            if kk_subj > 16:  # to avoid pb with subj 16 without shared data
                subject_id_bids = "{:02d}".format(
                    kk_subj - 1
                )  # zero padding to account for >10 subjects in this dataset
            else:
                subject_id_bids = "{:02d}".format(kk_subj)
            session_id = "{:02d}".format(kk_session)
            raw_eeg.pick(picks=["eeg"], exclude=["emg", "eog", "meg", "ecg"])
            print(np.unique(raw_eeg.get_channel_types()))  # to check the selection
            raw_eeg_list.append(raw)
            bids_path_eeg = BIDSPath(
                subject=subject_id_bids,
                session=session_id,
                task=task,
                run=run_id,
                root=bids_root,
                datatype="eeg",
            )
            bids_eeg_list.append(bids_path_eeg)
            write_raw_bids(
                raw=raw_eeg,
                bids_path=bids_path_eeg,
                events=events_array,
                event_id=event_dict,
                event_metadata=metadata,
                extra_columns_descriptions=metadata_description,
                montage=mon,
                overwrite=True,
                format="auto",
            )  # fif format does not work here...

            raw_meg.pick(picks=["meg"], exclude=["emg", "eog", "eeg", "ecg"])
            print(np.unique(raw_meg.get_channel_types()))  # to check the selection
            raw_meg_list.append(raw_meg)
            bids_path_meg = BIDSPath(
                subject=subject_id_bids,
                session=session_id,
                task=task,
                run=run_id,
                root=bids_root,
                datatype="meg",
            )
            bids_meg_list.append(bids_path_meg)
            write_raw_bids(
                raw=raw_meg,
                bids_path=bids_path_meg,
                events=events_array,
                event_id=event_dict,
                event_metadata=metadata,
                extra_columns_descriptions=metadata_description,
                overwrite=True,
                format="auto",
            )

counts_eeg = count_events(bids_root, datatype="eeg")
counts_eeg
counts_meg = count_events(bids_root, datatype="meg")
counts_meg
dataset_report = make_report(root=bids_root)
files2debug.to_csv(bids_root + "FileDebug_bidsify.csv")
