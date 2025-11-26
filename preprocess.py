import random
import sys
from Bio.Affy import CelFile

import pydicom
from pydicom import dcmread
from pydicom.data import get_testdata_file
import nibabel as nib

import numpy as np
import pandas as pd
import os
import gzip
import shutil
import zlib
import pysam
import matplotlib.pyplot as plt
from skimage.transform import resize

def unzip_gz(fp):
    with gzip.open(fp, 'rb') as f_in:
        with open(fp[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            os.remove(fp)
####################################### general clinical data ##########################################################################
# TODO: determine and encode categorical data
# TODO: standardize and normalize
########################################################################################################################################
def mimic_remove_males():
    patients = pd.read_csv(f'{path}/PATIENTS.csv')
    males = patients["SUBJECT_ID"][patients["GENDER"] == "M"].tolist() # contains SUBJECT_ID of male patients
    check_males = patients["SUBJECT_ID"][patients.GENDER == 'M'].tolist()
    print(len(check_males) == len(males)) # sanity check

    for file in os.listdir(path):
        fp = f'{path}/{file}'
        if not os.path.isdir(fp) and file.split(".")[-1] == "gz":
            with gzip.open(fp, 'rb') as f_in:
                with open(fp[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    os.remove(fp)
        elif file.split(".")[-1] == "csv" and file[-6:] != "_0.csv": # should not repeat infinitely
            cols = pd.read_csv(fp, nrows=0).columns.values.tolist() # pd.read_csv().columns.values returns np.ndarray
            new_fp = f'{path}/{file.split(".")[0]}_0.csv'
            if "SUBJECT_ID" in cols:
                df = pd.DataFrame(columns=cols)
                df.to_csv(new_fp, index=False, mode='w')
                for chunk in pd.read_csv(fp, chunksize=1000): 
                    print("chunking")
                    filtered = chunk[~chunk["SUBJECT_ID"].isin(males)] # only reads subject id's corresponding to female patients
                    filtered.to_csv(new_fp, index=False, mode='a', header=False)
            # csv's without subject information: D_CPT, D_ITEMS, CAREGIVERS, D_LABITEMS, D_ICD_DIAGNOSES, D_ICD_PROCEDURES
        else:
            print(file)

def mimic_drop_duplicates():
    # drop duplicates
    for file in os.listdir(path):
        if file[-6:] == "_0.csv":
            fp = f'{path}/{file}'
            print(file)
            df = pd.read_csv(fp, low_memory=False) # pyarrow reads csv rows poorly and mistakes commas within entries as a new column
            initial = len(df["SUBJECT_ID"])
            df = df.drop_duplicates()
            if len(df["SUBJECT_ID"]) < initial:
                print("writing")
                df.to_csv(fp)
            elif len(df["SUBJECT_ID"]) == initial:
                print("do nothing")
            else:
                print("something is wrong")

path = "/Users/vanigupta/Documents/uf_digital_twin/mimic-iii-clinical-database-1.4"
# mimic_remove_males()
# mimic_drop_duplicates()
####################################### hormonal time series data ######################################################################
# made my own data :p
########################################################################################################################################
def generate_hormonal_timeseries():
    np.random.seed(42)

    n_users = 100

    user_ids = [f"user_{i+1}" for i in range(n_users)]
    cycle_lengths = np.random.randint(26, 37, size=n_users)

    diagnoses = np.random.choice(["PCOS", "fibroids", "normal"], size=n_users, p=[0.2, 0.3, 0.5])

    records = []
    for user_id, cycle_len, dx in zip(user_ids, cycle_lengths, diagnoses):
        menstrual_len = np.random.randint(4,7)
        follicular_len = menstrual_len + np.random.randint(7,14)
        ovulation_len = follicular_len + np.random.randint(3,5)
        remaining_luteal = cycle_len - ovulation_len
        early_luteal_len = ovulation_len + int(0.5 * remaining_luteal)
        late_luteal_len = cycle_len

        for day in range(1, cycle_len + 1):
            # simulate hormone values based on phase
            if day <= menstrual_len:
                estradiol = 20 + np.random.normal(0, 2)
                progesterone = 0.3 + np.random.normal(0, 0.1)
                phase = "menstruation"
            elif day <= follicular_len:
                estradiol = 30 + 10 * np.sin(day / follicular_len * np.pi) + np.random.normal(0, 3)
                progesterone = 0.5 + np.random.normal(0, 0.2)
                phase = "follicular"
            elif day <= ovulation_len:
                phase = "ovulation"
                estradiol = 100 + np.random.normal(0, 5)
                progesterone = 1.5 + np.random.normal(0, 0.5)
            elif day <= early_luteal_len:
                phase = "early_luteal"
                estradiol = 60 + np.random.normal(0, 3)
                progesterone = 10 + np.random.normal(0, 1)
            else:
                phase = "late_luteal"
                estradiol = 40 + np.random.normal(0, 2)
                progesterone = 6 + np.random.normal(0, 1)

            # fibroid-specific elevated estrogen plateau
            if dx == "fibroids" and phase in ["follicular", "menstruation"]:
                estradiol += 10

            estrone = 0.5 * estradiol + np.random.normal(0, 2)
            testosterone = 0.4 + 0.1 * np.sin((day - 7) / cycle_len * 2 * np.pi) + np.random.normal(0, 0.02)
            hcg = 0 + np.random.normal(0, 0.05)

            records.append({
                "User ID": user_id,
                "Day": day,
                "Cycle Length": cycle_len,
                "Cycle Phase": phase,
                "Diagnosis": dx,
                "Estradiol (E2)": estradiol,
                "Estrone (E1)": estrone,
                "Progesterone": progesterone,
                "Testosterone": testosterone,
                "HCG": hcg
            })

    df = pd.DataFrame(records)
    output_path = "simulated_hormone_cycles.csv"
    df.to_csv(output_path, index=False)
####################################### uf-specific data ###############################################################################
# genetic data
########################################################################################################################################
# path = "/Users/vanigupta/Documents/uf_digital_twin/GSE30673_RAW"
# for file in os.listdir(path):
#     fp = f'{path}/{file}'
#     if not os.path.isdir(fp) and file.split(".")[-1] == "gz":
#         unzip_gz(fp)
#     elif file.split(".")[-1] == "CEL":
#         try:
#             with open(fp, 'rb') as cel_file:
#                 c = CelFile.read(cel_file)
#         except:
#             print(f"{file} not working") # GSM760701.CEL, GSM760703.CEL
#     else:
#         print(file)

# ####################################### mri imaging data ###############################################################################
# # DICOM header file contains this info: (a) Patient (b) Study (c) Series (d) Image
# # seg files contain the label: (1) uterine wall, (2) uterine cavity, (3) myoma, or (4) nabothian cyst
# ########################################################################################################################################
def extract_mri_data():
    path = "/Users/vanigupta/Documents/uf_digital_twin/UMD"
    count = 0
    patient_records = []

    for patient_id in os.listdir(path):
        patient_path = os.path.join(path, patient_id)
        if not os.path.isdir(patient_path) or patient_id.startswith('.'):
            continue
        
        patient_id = patient_id.split("_")[-1]
        record = {
            "patient_id": patient_id,
            "patient_weight": None,
            "labels": [],
            "fibroid_present": False,
            "downsampled_shape": None
        }

        for file in os.listdir(patient_path):
            fp = os.path.join(patient_path, file)

            if file.endswith(".dcm") and record["patient_weight"] is None:
                try:
                    ds = pydicom.dcmread(fp, stop_before_pixels=True)
                    record["patient_weight"] = float(getattr(ds, "PatientWeight", None))
                except:
                    pass
            elif file.endswith("_t2.nii"):
                try:
                    t2 = nib.load(fp)
                    t2_data = t2.get_fdata()
                except:
                    pass
            elif file.endswith("_seg.nii"):
                try:
                    seg = nib.load(fp)
                    seg_data = seg.get_fdata()
                    labels = np.unique(seg_data).astype(int).tolist()
                    record["labels"] = labels
                    if 3 in labels:
                        record["fibroid_present"] = True
                except:
                    pass

        # downsample
        if t2_data is not None and seg_data is not None:
            try:
                target_shape = (128, 128, 5)
                t2_down = resize(t2_data, target_shape, order=1, anti_aliasing=True)
                seg_down = resize(seg_data, target_shape, order=0, preserve_range=True).astype(np.uint8)

                # z-score normalization
                std = np.std(t2_down)
                # avoid division by zero in case the volume is constant
                if std > 1e-6:
                    t2_down = (t2_down - np.mean(t2_down)) / std
                else:
                    t2_down = t2_down - np.mean(t2_down)  # mean-center if std is too small

                np.save(f"{patient_id}_t2_downsampled.npy", t2_down)
                np.save(f"{patient_id}_seg_downsampled.npy", seg_down)

                record["downsampled_shape"] = list(t2_down.shape)
            except Exception as e:
                print(f"[DOWNSAMPLE ERROR] {patient_id}: {e}")

        patient_records.append(record)
        count += 1
        print(count)

    df = pd.DataFrame(patient_records)
    df.to_csv("umd_data_output.csv", index=False)


generate_hormonal_timeseries()
extract_mri_data()
# df_small = pd.read_csv("mri_vectors.csv", nrows=1000)
# df_small.to_csv("mri_vectors_small.csv", index=False)