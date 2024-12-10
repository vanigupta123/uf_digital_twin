import pydicom
from pydicom import dcmread
from pydicom.data import get_testdata_file
import numpy as np
import nibabel as nib
import os
from sh import gunzip

####################################### general clinical data ##########################################################################
# TODO: remove all male patients
# TODO: handle any missing or duplicate data
# TODO: determine and encode categorical data
# TODO: standardize and normalize
########################################################################################################################################

####################################### hormonal time series data ######################################################################
# TODO: handle any missing or duplicate data
# TODO: determine and encode categorical data
# TODO: standardize and normalize
# TODO: aggregate data and resample time intervals
########################################################################################################################################

####################################### uf-specific data ###############################################################################

########################################################################################################################################

####################################### mri imaging data ###############################################################################
# DICOM header file contains this info: (a) Patient (b) Study (c) Series (d) Image
# seg files contain the label: (1) uterine wall, (2) uterine cavity, (3) myoma, or (4) nabothian cyst
########################################################################################################################################
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
                gunzip(fp, f"-c > {fp[:-3]}")
        elif file.split(".")[-1] == "nii":
            img = nib.load(fp).get_fdata()
        else:
            print(f'{fp} unhandled')
