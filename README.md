# Fetal-brain localization and orientation estimation

![Fetal-brain geometry reconstruction](fetal-align.png)

## Purpose

This project comprises fully automated brain localization, extraction and
orientation estimation for fetal-brain MRI. The estimation is achieved in
seconds using blob-detection techniques to find the brain and eyes. 

## Prerequisites

MATLAB version 9.1/R2016b or later is required. If necessary, the pre-compiled
MEX function can be rebuilt by running `mex mser.cpp` in MATLAB.

## Test data

The test data includes 41 EPI and 2 HASTE acquisitions. Except for
[ep2d_34.nii.gz](data/ep2d_34.nii.gz) and [ep2d_41.nii.gz](data/ep2d_34.nii.gz),
the fetal brain and eyes are correctly detected by the algorithm.

## Example

For a demo showcasing the different stages of the algorithm, run `demo`. The
`printstats` script compares the landmarks derived by the algorithm to manually
localized eye and brain centers.

## FreeSurfer

The IO scripts included in [freesurfer](freesurfer) are distributed under the
FreeSurfer Software License Agreement. See https://surfer.nmr.mgh.harvard.edu
and [freesurfer/LICENSE](freesurfer/LICENSE) for details.

## References

Manuscript under review.
