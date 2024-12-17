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

####################################### general clinical data ##########################################################################
# TODO: determine and encode categorical data
# TODO: standardize and normalize
########################################################################################################################################
# path = "/Users/vanigupta/Documents/uf_digital_twin/mimic-iii-clinical-database-1.4"
# patients = pd.read_csv(f'{path}/PATIENTS.csv')
# males = patients["SUBJECT_ID"][patients["GENDER"] == "M"].tolist() # contains SUBJECT_ID of male patients
# check_males = patients["SUBJECT_ID"][patients.GENDER == 'M'].tolist()
# print(len(check_males) == len(males)) # sanity check

# for file in os.listdir(path):
#     fp = f'{path}/{file}'
#     if not os.path.isdir(fp) and file.split(".")[-1] == "gz":
#         with gzip.open(fp, 'rb') as f_in:
#             with open(fp[:-3], 'wb') as f_out:
#                 shutil.copyfileobj(f_in, f_out)
#                 os.remove(fp)
#     elif file.split(".")[-1] == "csv" and file[-6:] != "_0.csv": # should not repeat infinitely
#         cols = pd.read_csv(fp, nrows=0).columns.values.tolist() # pd.read_csv().columns.values returns np.ndarray
#         new_fp = f'{path}/{file.split(".")[0]}_0.csv'
#         if "SUBJECT_ID" in cols:
#             df = pd.DataFrame(columns=cols)
#             df.to_csv(new_fp, index=False, mode='w')
#             for chunk in pd.read_csv(fp, chunksize=1000): 
#                 print("chunking")
#                 filtered = chunk[~chunk["SUBJECT_ID"].isin(males)] # only reads subject id's corresponding to female patients
#                 filtered.to_csv(new_fp, index=False, mode='a', header=False)
#         # csv's without subject information: D_CPT, D_ITEMS, CAREGIVERS, D_LABITEMS, D_ICD_DIAGNOSES, D_ICD_PROCEDURES
#     else:
#         print(file)

# # drop duplicates
# for file in os.listdir(path):
#     if file[-6:] == "_0.csv":
#         fp = f'{path}/{file}'
#         print(file)
#         df = pd.read_csv(fp, low_memory=False) # pyarrow reads csv rows poorly and mistakes commas within entries as a new column
#         initial = len(df["SUBJECT_ID"])
#         df = df.drop_duplicates()
#         if len(df["SUBJECT_ID"]) < initial:
#             print("writing")
#             df.to_csv(fp)
#         elif len(df["SUBJECT_ID"]) == initial:
#             print("do nothing")
#         else:
#             print("something is wrong")
####################################### hormonal time series data ######################################################################
# TODO: handle any missing or duplicate data
# TODO: determine and encode categorical data
# TODO: standardize and normalize
# TODO: aggregate data and resample time intervals
########################################################################################################################################

####################################### uf-specific data ###############################################################################
# genetic data
########################################################################################################################################
# path = "/Users/vanigupta/Documents/uf_digital_twin/GSE30673_RAW"
# for file in os.listdir(path):
#     fp = f'{path}/{file}'
#     if not os.path.isdir(fp) and file.split(".")[-1] == "gz":
#         with gzip.open(fp, 'rb') as f_in:
#             with open(fp[:-3], 'wb') as f_out:
#                 shutil.copyfileobj(f_in, f_out)
#                 os.remove(fp)
#     elif file.split(".")[-1] == "CEL":
#         with open(fp, 'rb') as cel_file:
#             c = CelFile.read(cel_file) # debug this next
#     else:
#         print(file)

# ####################################### mri imaging data ###############################################################################
# # DICOM header file contains this info: (a) Patient (b) Study (c) Series (d) Image
# # seg files contain the label: (1) uterine wall, (2) uterine cavity, (3) myoma, or (4) nabothian cyst
# ########################################################################################################################################
path = "/Users/vanigupta/Documents/uf_digital_twin/UMD"
# delete = list()
delete = ['113_seq.nii', '050_seq.nii', '061_seq.nii', '059_seq.nii', '200_seq.nii', '060_seq.nii', '241_seq.nii', 
'074_seq.nii', '045_seq.nii', '222_seq.nii', '249_seq.nii', '088_seq.nii', '086_seq.nii', '130_seq.nii', '106_seq.nii', 
'139_seq.nii', '164_seq.nii', '127_seq.nii', '054_seq.nii', '062_seq.nii', '037_seg.nii', '063_seq.nii', '064_seq.nii', 
'090_seq.nii', '055_seq.nii', '070_seq.nii', '079_seq.nii', '046_seq.nii', '085_seq.nii', '071_seq.nii', '076_seq.nii', 
'195_seq.nii', '157_seq.nii', '156_seq.nii']
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
        # elif file.split(".")[-1] == "gz": # no need to run/unzip more than once
        #     try:
        #         with gzip.open(fp, 'rb') as f_in:
        #             with open(fp[:-3], 'wb') as f_out:
        #                 shutil.copyfileobj(f_in, f_out)
        #     except Exception as error:
        #         with open(fp, "rb") as f:
        #             content = f.read(16)
        #             if '000000000' not in content.hex(): # check if .gz file is salvagable
        #                 print(file)
        #             else:
        #                 num = file.split("_")[-2]
        #                 idx = file.index(num)
        #                 delete.append(file[idx:-3])
        elif not os.path.isdir(fp) and file.split(".")[-1] == "nii":
            # print(file[-11:])
            if file[-11:] in delete: # make sure there are no longer any empty .nii files
                print(f'deleting {file}')
                os.remove(fp)
                continue
            try:
                nii = nib.load(fp)
                img = nii.get_fdata()
            except Exception as error: # find remaining empty fles or other exceptions
                print(f'{file}: {error}')
            if file[-8:] == "_seg.nii": # retrieve labels
                # print(nii.header)
                print(np.unique(img))
        # else: # just folders, .gz, and .DStore
        #     print(f'{fp} unhandled')