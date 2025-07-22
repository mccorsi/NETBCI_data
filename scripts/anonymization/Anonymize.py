##% Code to anonymize the BIDS dataset
# Author: Marie-Constance Corsi
# Last modification: July, 22nd 2025

# %%
import os
import os.path as op
import shutil

import pandas as pd
from mne_bids import (
    anonymize_dataset,
    print_dir_tree,
)

# %% Paths
if os.path.basename(os.getcwd()) == "scripts":
    os.chdir("..")
basedir = os.getcwd()

os.chdir(basedir)
bids_root = basedir + "/NETBCI_postBIDS_db/"
warp_path = op.join(bids_root, "derivatives", "warpimg")
bids_root_anon = bids_root + "/NETBCI_postBIDS_db_anon/"
path_figures_root = basedir + "/Figures/"

# %% Anonymize the dataset
anonymize_dataset(
    bids_root_in=bids_root,
    bids_root_out=bids_root_anon,
    datatypes=["meg", "eeg", "anat"],
)

print_dir_tree(bids_root_anon)

# %% Save somewhere the correspondence between pre-post shuffle of the subjects' ID:
dict_ids = dict()
list_shuffle = (
    []
)  # list of the subjects_ID post shuffle displayed as output of the previous section
anonymized_id = ["sub-{:02d}".format(kk_subj) for kk_subj in list_shuffle]
subject_id = ["sub-{:02d}".format(kk_subj) for kk_subj in range(1, 20)]
dict_ids["subject_id"] = subject_id
dict_ids["anonymized_id"] = anonymized_id
pd_ids = pd.DataFrame.from_dict(dict_ids)
pd_ids.to_csv(bids_root + "linkdata.csv")

# %% Update participants.csv file accordingly
df_participants = pd.read_csv(bids_root_anon + "participants.tsv", sep="\t")
tmp_df_participants = df_participants["participant_id"]

df_participants_orig = pd.read_csv(bids_root + "participants.tsv", sep="\t")
df_participants_orig_copy = df_participants_orig.copy()
df_participants_orig_copy["Unnamed: 0"] = [
    "sub-{:02d}".format(kk_subj) for kk_subj in list_shuffle
]

df_participants_orig_copy = df_participants_orig_copy.sort_values(by=["Unnamed: 0"])
df_participants_orig_copy2 = df_participants_orig_copy.copy()
df_participants_orig_copy2.drop(["participant_id"], axis=1, inplace=True)
df_participants_annon = df_participants_orig_copy2
df_participants_annon = df_participants_annon.rename(
    columns={"Unnamed: 0": "participant_id"}
)

df_participants_annon.to_csv(bids_root_anon + "participants.tsv", sep="\t")


# %% Rearrange warped MRI accordingly - adapted from the code proposed by @mikkel
overwrite = 0
nomri_org = []
nomri_ano = []

# Make derivatives directory
drv_path = op.join(bids_root_anon, "derivatives", "warpimg")
if not op.exists(drv_path):
    os.makedirs(drv_path)

# Copy data
linkdata = pd.read_csv(bids_root + "linkdata.csv")
for ii, subjid in enumerate(linkdata.anonymized_id):
    origid = linkdata.subject_id[ii]

    for kk_sess in range(1, 5):
        sess_id = "ses-{:02d}".format(kk_sess)

        infile = op.join(
            warp_path,
            origid,
            sess_id,
            "anat",
            origid + "_" + sess_id + "_mri_warptmp.nii",
        )

        out_path = op.join(drv_path, subjid, sess_id, "anat")
        outfile = op.join(out_path, subjid + "_" + sess_id + "_desc-warpimg_T1w.nii")

        if not op.exists(infile):
            nomri_org += [origid]
            nomri_ano += [subjid]
            continue

        if op.exists(outfile) and not overwrite:
            print("Output exist for subj " + subjid)
            continue

        if not op.exists(out_path):
            os.makedirs(out_path)

        shutil.copy2(infile, outfile)
