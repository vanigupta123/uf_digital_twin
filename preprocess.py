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

####################################### general clinical data ##########################################################################
# TODO: determine and encode categorical data
# TODO: standardize and normalize
########################################################################################################################################
path = "/Users/vanigupta/Documents/uf_digital_twin/mimic-iii-clinical-database-1.4"
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
####################################### hormonal time series data ######################################################################
# TODO: handle any missing or duplicate data
# TODO: determine and encode categorical data
# TODO: standardize and normalize
# TODO: aggregate data and resample time intervals
########################################################################################################################################

####################################### uf-specific data ###############################################################################
# genetic data
########################################################################################################################################
path = "/Users/vanigupta/Documents/uf_digital_twin/GSE30673_RAW"
for file in os.listdir(path):
    fp = f'{path}/{file}'
    if not os.path.isdir(fp) and file.split(".")[-1] == "gz":
        with gzip.open(fp, 'rb') as f_in:
            with open(fp[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                os.remove(fp)
    elif file.split(".")[-1] == "CEL":
        with open(fp, 'rb') as cel_file:
            c = CelFile.read(cel_file, version=4)
    else:
        print(file)
    

# ####################################### mri imaging data ###############################################################################
# # DICOM header file contains this info: (a) Patient (b) Study (c) Series (d) Image
# # seg files contain the label: (1) uterine wall, (2) uterine cavity, (3) myoma, or (4) nabothian cyst
# ########################################################################################################################################
path = "/Users/vanigupta/Documents/uf_digital_twin/UMD"
for folder in os.listdir(path):
    next_path = f'{path}/{folder}'
    if folder == ".DS_Store":
        continue
    for file in os.listdir(next_path):
        fp = os.path.join(next_path, file)
        if file == ".DS_Store":
            continue
        elif file.split(".")[-1] == "dcm":
            ds = dcmread(fp, force=True)
            # may need to read header
            # header = pydicom.filereader.read_preamble(fp, force=True)
            # print(header)
        elif file.split(".")[-1] == "gz": # unzips nii.gz files
            with gzip.open(fp, 'rb') as f_in:
                with open(fp[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        elif not os.path.isdir(fp) and file.split(".")[-1] == "nii":
            img = nib.load(fp).get_fdata()
        else:
            print(f'{fp} unhandled')
