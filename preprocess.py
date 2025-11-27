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
from sklearn.preprocessing import StandardScaler
from skimage.transform import resize
from skimage.measure import label

def unzip_gz(fp):
    with gzip.open(fp, 'rb') as f_in:
        with open(fp[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            os.remove(fp)
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
                "user id": user_id,
                "day": day,
                "cycle length": cycle_len,
                "cycle phase": phase,
                "diagnosis": dx,
                "estradiol (E2)": estradiol,
                "estrone (E1)": estrone,
                "progesterone": progesterone,
                "testosterone": testosterone,
                "HCG": hcg
            })

    df = pd.DataFrame(records)

    # normalize
    hormone_cols = ["estradiol (E2)", "estrone (E1)", "progesterone", "testosterone", "HCG"]
    df[hormone_cols] = (df[hormone_cols] - df[hormone_cols].mean()) / df[hormone_cols].std()

    df_encoded = pd.get_dummies(df, columns=["cycle phase", "diagnosis"], dummy_na=True)


    output_path = "simulated_hormone_cycles.csv"
    df_encoded.to_csv(output_path, index=False)
# ####################################### mri imaging data ###############################################################################
# # DICOM header file contains this info: (a) Patient (b) Study (c) Series (d) Image
# # seg files contain the label: (1) uterine wall, (2) uterine cavity, (3) myoma, or (4) nabothian cyst
# ########################################################################################################################################
def generate_synthetic_metadata(patients):
    synthetic_records = []
    age_order = {"18–29": 0, "30–44": 1, "45+": 2}
    for patient in patients:
        fib_count = patient.get("num_fibroids", 0)
        fib_ratio = patient.get("fibroid_volume_ratio", 0.0)

        if fib_count == 0:
            pain = 0
        elif fib_ratio < 0.4:
            pain = min(8, int(np.random.rand()*fib_ratio*100))
        else:
            pain = random.randint(7, 10)

        if fib_count == 0:
            treatment = "none"
        elif fib_ratio >= 0.4:
            treatment = "surgery"
        elif (fib_ratio > 0.2 and fib_ratio < 0.4) or pain > 5:
            treatment = random.choices(["surgery", "hormonal"], weights=[0.7, 0.3])[0]
        else:
            treatment = random.choices(["hormonal", "none"], weights=[0.8, 0.2])[0]

        # if hormonal
        hormonal_mods = ["estrogen decreased", "progesterone increased", "GnRH agonist"]
        hormonal_mod = random.choice(hormonal_mods) if treatment == "hormonal" else None

        # is this relevant?
        ethnicity = random.choices(
            ["White", "Black", "Asian", "Hispanic", "Other"],
            weights=[0.25, 0.25, 0.20, 0.20, 0.10]
        )[0]

        age_group = random.choices(["18–29", "30–44", "45+"], weights=[0.3, 0.5, 0.2])[0]

        prior_preg = random.random() < 0.6 if age_group != "18–29" else random.random() < 0.3

        patient.update({
            "pain_level": pain,
            "treatment_type": treatment,
            "hormonal_mod": hormonal_mod,
            "ethnicity": ethnicity,
            "age_group": age_group,
            "age_group_encoded": age_order[age_group],
            "prior_pregnancy": prior_preg
        })
        synthetic_records.append(patient)

    df = pd.DataFrame(synthetic_records)

    # normalize
    scaler = StandardScaler()
    df[["pain_level", "num_fibroids", "fibroid_volume_ratio"]] = scaler.fit_transform(
        df[["pain_level", "num_fibroids", "fibroid_volume_ratio"]]
    )

    df_encoded = pd.get_dummies(df, columns=[
        "treatment_type", "ethnicity", "prior_pregnancy", "hormonal_mod"
    ], dummy_na=True)

    return df_encoded

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
            "num_fibroids": 0.0,
            "fibroid_volume_ratio": 0.0,
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
                    
                    fibroid_voxels = (seg_data == 3)
                    labeled_regions = label(fibroid_voxels)
                    num_fibroids = labeled_regions.max()  # number of disconnected regions / clusters

                    fibroid_volume = np.sum(seg_data == 3)
                    other_volume = np.sum((seg_data > 0) & (seg_data != 3))
                    ratio = fibroid_volume / (fibroid_volume + other_volume + 1e-6)

                    record["num_fibroids"] = num_fibroids
                    record["fibroid_volume_ratio"] = ratio

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

                np.save(f"umd_np_arrays/{patient_id}_t2_downsampled.npy", t2_down)
                np.save(f"umd_np_arrays/{patient_id}_seg_downsampled.npy", seg_down)

                record["downsampled_shape"] = list(t2_down.shape)
            except Exception as e:
                print(f"[DOWNSAMPLE ERROR] {patient_id}: {e}")

        patient_records.append(record)
        count += 1

    df = generate_synthetic_metadata(patient_records)
    df.to_csv("umd_data_categorical.csv", index=False)

# ########################################################################################################################################

generate_hormonal_timeseries()
extract_mri_data()