##% Code to update metadata related to the participants
# Author: Marie-Constance Corsi
# Last modification: July, 22nd 2025

# %%
import os

import pandas as pd

# %%
if os.path.basename(os.getcwd()) == "scripts":
    os.chdir("..")
basedir = os.getcwd()

os.chdir(basedir)
path_data2bidsify_root = basedir + "/NETBCI_db/"
path_dataPostbidsify_root = basedir + "/NETBCI_postBIDS_db/"
list_mod = ["MEG", "EEG"]

bids_root = path_dataPostbidsify_root

# for each subject and each session list sub-0X_ses-0X_scans.tsv files & merge infos
list_subj = list(range(1, 20))

# %% Put the list of the generated files in a single doc to be used afterwards
tsv = pd.DataFrame()

for kk_subj in list_subj:
    subject_id = "{:02d}".format(
        kk_subj
    )  # zero padding to account for >10 subjects in this dataset

    for kk_session in list(range(1, 5)):
        session_id = "{:02d}".format(kk_session)
        path_subj_sess = (
            path_dataPostbidsify_root + "sub-" + subject_id + "/ses-" + session_id
        )

        for file in os.listdir(path_subj_sess):

            if file.endswith("scans.tsv"):
                tmp = pd.read_csv(os.path.join(path_subj_sess, file), sep="\t")
                tsv = pd.concat([tsv, tmp], ignore_index=True)
                tsv.to_csv(
                    path_dataPostbidsify_root + "/filenames.csv",
                    sep=",",
                    header=True,
                    # sep="\t", header=True,
                    index=False,
                )

# %% Dataset description update -> cf create_metadata_filenames.m & update_infos4BIDS.m

# %% update participants description json file
import json

with open("replayScript.json", "r+") as jsonFile:
    data = json.load(jsonFile)

    data["location"] = "NewPath"

    jsonFile.seek(0)  # rewind
    json.dump(data, jsonFile)
    jsonFile.truncate()


# %% update Participants.tsv with all the information even for each session - information inserted manually before the anonymization
participantsV1 = pd.read_csv(path_dataPostbidsify_root + "participants.tsv", sep="\t")
participantsV1.drop(["height", "weight"], axis=1, inplace=True)

participants = pd.DataFrame()
participants["participant_id"] = [
    "sub-" + "{:02d}".format(kk_subj) for kk_subj in range(1, 20)
]
participants["age"] = []
participants["sex"] = []
participants["hand"] = ["R"] * len(participants)
participants["RosenbergScale"] = []
participants["EMG28-amotivation"] = []
participants["EMG28-extrinsic-externReg"] = []
participants["EMG28-extrinsic-identifReg"] = []
participants["EMG28-extrinsic-introjReg"] = []
participants["EMG28-intrinsic.accomplishment"] = []
participants["EMG28-intrinsic.knowledge"] = []
participants["EMG28-intrinsic.stimulation"] = []
participants["VMIQ2-externalVisualImagery"] = []
participants["VMIQ2-InternalVisualImagery"] = []
participants["VMIQ2-KinaestheticImagery"] = []
participants["STAI_YA-session1"] = []
participants["STAI_YA-session2"] = []
participants["STAI_YA-session3"] = []
participants["STAI_YA-session4"] = []

participants.to_csv(path_dataPostbidsify_root + "participants.tsv", sep="\t")
