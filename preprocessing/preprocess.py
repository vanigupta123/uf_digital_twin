import random
import sys
# from Bio.Affy import CelFile

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
    # hormone_cols = ["estradiol (E2)", "estrone (E1)", "progesterone", "testosterone", "HCG"]
    # df[hormone_cols] = (df[hormone_cols] - df[hormone_cols].mean()) / df[hormone_cols].std()
    # shouldn't normalize before splitting -- will cause data leakage / distribution shift

    df_encoded = pd.get_dummies(df, columns=["cycle phase", "diagnosis"], dummy_na=True)


    output_path = "simulated_hormone_cycles.csv"
    df_encoded.to_csv(output_path, index=False)

# ####################################### generate synthetic metadata #######################################################################
def generate_synthetic_metadata(patients):
    synthetic_records = []
    age_order = {"18–29": 0, "30–44": 1, "45+": 2}
    for patient in patients:
        fib_count = patient.get("num_fibroids", 0)
        fib_ratio = patient.get("fibroid_volume_ratio", 0.0)

        if fib_count == 0:
            pain = 0
        elif fib_ratio < 0.4:
            # mild-moderate: pain skewed toward lower values but occasionally higher
            pain = int(np.random.beta(2, 5) * 10)
        else:
            # severe: pain skewed toward higher values
            pain = int(np.random.beta(5, 2) * 10)
        pain = min(10, max(0, pain))

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

        # black women more likely to have fibroids
        ethnicity = random.choices(
            ["White", "Black", "Asian", "Hispanic", "Other"],
            weights=[0.25, 0.25, 0.20, 0.20, 0.10]
        )[0]

        age_group = random.choices(["18–29", "30–44", "45+"], weights=[0.3, 0.5, 0.2])[0]

        prior_preg = random.random() < 0.6 if age_group != "18–29" else random.random() < 0.3

        # ferritin_proxy: missing ~40% of rows (requires blood draw)
        if random.random() < 0.4:
            patient["ferritin_proxy"] = None  # genuinely missing / no blood draw
        else:
            # low ferritin is associated with heavy bleeding / fibroids
            patient["ferritin_proxy"] = round(random.gauss(
                25 if fib_ratio > 0.3 else 45, 12
            ), 1)

        if random.random() < 0.05:
            patient["cycle_length_days"] = None
        elif fib_ratio > 0.3:
            patient["cycle_length_days"] = round(random.gauss(35, 5), 0)
        else:
            patient["cycle_length_days"] = round(random.gauss(28, 5), 0)
        patient["cycle_length_days"] = None if patient["cycle_length_days"] is None else min(60, max(18, patient["cycle_length_days"]))

        if random.random() < 0.25:
            patient["flow_intensity"] = None
        elif fib_ratio > 0.3:
            patient["flow_intensity"] = random.choices(["moderate", "heavy", "very_heavy"], weights=[0.1, 0.5, 0.4])[0]
        else:
            patient["flow_intensity"] = random.choices(["light", "moderate", "heavy"], weights=[0.2, 0.5, 0.3])[0]
        
        if random.random() < 0.3:
            patient["symptom_duration_months"] = None
        elif patient["fibroid_present"]:
            patient["symptom_duration_months"] = random.gauss(18, 12)
            patient["symptom_duration_months"] = min(120, max(0, patient["symptom_duration_months"]))
        else:
            patient["symptom_duration_months"] = 0

        patient.update({
            "pain_level": pain,
            "treatment_type": treatment,
            # "hormonal_mod": hormonal_mod,
            "ethnicity": ethnicity,
            "age_group": age_group,
            "age_group_encoded": age_order[age_group],
            "prior_pregnancy": prior_preg
        })
        synthetic_records.append(patient)

    df = pd.DataFrame(synthetic_records)
    df["flow_intensity"] = df["flow_intensity"].map({"light": 1, "moderate": 2, "heavy": 3, "very_heavy": 4})
    df["flow_intensity"] = df["flow_intensity"].fillna(0)
    df["flow_intensity_missing"] = df["flow_intensity"] == 0
    df_encoded = pd.get_dummies(df, columns=[
        "treatment_type", "ethnicity", 
    # "prior_pregnancy", "hormonal_mod", 
    #     "flow_intensity"
    ])

    return df_encoded
# ####################################### mri imaging data ###############################################################################
# # DICOM header file contains this info: (a) Patient (b) Study (c) Series (d) Image
# # seg files contain the label: (1) uterine wall, (2) uterine cavity, (3) myoma, or (4) nabothian cyst
# ########################################################################################################################################

def extract_mri_data():
    path = "/Users/vanigupta/Documents/uf_digital_twin/UMD"
    count = 0
    patient_records = []
    t2_data = None
    seg_data = None

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
                except Exception as e:
                    pass
            elif file.endswith("_t2.nii"):
                try:
                    t2 = nib.load(fp)
                    t2_data = t2.get_fdata()
                except Exception as e:
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

                except Exception as e:
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
def fibroid_growth_dataset():
    cols = ["patient_id","patient_weight","num_fibroids","fibroid_volume_ratio","flow_intensity","symptom_duration_months","pain_level","age_group","age_group_encoded","prior_pregnancy","flow_intensity_missing","treatment_type_hormonal","treatment_type_none","treatment_type_surgery","ethnicity_Asian","ethnicity_Black","ethnicity_Hispanic","ethnicity_Other","ethnicity_White"]
    # add time interval column and maybe join this dataset with the categorical dataset by patient id
    # example dataset
    # [patient_id="92", 30 days, num_fibroids=5, fibroid_volume_ratio=0.76, patient_weight=60kg, flow_intensity=high, pain_level=4, treatment="blah"] 
    # maybe like 5-10 rows per patient with different time intervals. the time column can be in units of days
    # use the verhulst equation with the decay modifier to make up values for this i think
    # dV/dt = r * V * (1 - V/K) - d * V
    # closed form: V(t) = K / (1 + ((K - V0) / V0) * exp(-(r - d) * t))
    df = pd.read_csv("umd_data_categorical.csv")
    
    # growth/decay parameters
    r = 0.008  # base daily growth rate (slow — fibroids grow over months/years)
    K = 1.0    # carrying capacity (max volume ratio)
    
    # treatment-dependent decay rates
    # d > r means shrinkage, d < r means slower growth, d = 0 means no treatment effect
    decay_rates = {
        "none": 0.0,
        "hormonal": 0.012,   # slightly greater than r → slow shrinkage
        "surgery": 0.035     # much greater than r → fast shrinkage
    }
    
    # time points in days — unevenly spaced like real clinical checkups
    all_time_points = [0.0, 30.0, 60.0, 90.0, 120.0, 180.0, 270.0, 365.0, 540.0, 730.0]
    
    rows = []
    
    for _, patient in df.iterrows():
        pid = patient["patient_id"]
        V0 = patient["fibroid_volume_ratio"]
        
        # skip patients with no fibroids — nothing to model
        if V0 <= 0.001:
            continue
        
        # figure out treatment type from one-hot columns
        if patient.get("treatment_type_surgery", 0) == 1:
            treatment = "surgery"
        elif patient.get("treatment_type_hormonal", 0) == 1:
            treatment = "hormonal"
        else:
            treatment = "none"
        
        d = decay_rates[treatment]
        effective_rate = r - d  # positive = growth, negative = shrinkage
        
        # pick which time points this patient has observations at
        # randomly keep 3-5 time points per patient (simulates sparse clinical data)
        n_obs = random.randint(3, 5)
        time_points = sorted(random.sample(all_time_points, n_obs))
        # always include t=0 as the initial observation
        if 0 not in time_points:
            time_points = [0] + time_points[:n_obs - 1]
        
        for t in time_points:
            # modified Verhulst closed-form solution
            if abs(effective_rate) < 1e-10:
                # edge case: r ≈ d, no net growth or shrinkage
                V_t = V0
            else:
                denom = 1 + ((K - V0) / V0) * np.exp(-effective_rate * t)
                V_t = K / denom
            
            # add Gaussian noise to simulate measurement error (ultrasound isn't perfect)
            noise = np.random.normal(0, 0.02)
            V_t = np.clip(V_t + noise, 0.0, 1.0)
            
            rows.append({
                "patient_id": pid,
                "t_days": t,
                "fibroid_volume_ratio": round(V_t, 4),
                "num_fibroids": patient["num_fibroids"],
                "patient_weight": patient.get("patient_weight", None),
                "pain_level": float(patient.get("pain_level", None)),
                "age_group_encoded": patient.get("age_group_encoded", None),
                "treatment_type": treatment,
                "effective_rate": round(effective_rate, 4)
            })
    
    result = pd.DataFrame(rows)
    result.to_csv("pinn_fibroid_growth.csv", index=False)
    print(f"generated {len(result)} rows for {result['patient_id'].nunique()} patients")
    print(f"columns: {list(result.columns)}")
    print(f"sample:\n{result.head(10)}")
    return result

# extract_mri_data()
fibroid_growth_dataset()
# generate_hormonal_timeseries()