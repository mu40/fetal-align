// Extract maximally stable extremal regions (MSERs) from a 2D or 3D image (see
// Matas, J et al. Robust Wide Baseline Stereo from Maximally Stable Extremal
// Regions. 2002). Implementation inspired by Kristensen, F et al. Real-Time
// Extraction of Maximally Stable Extremal Regions on an FPGA. 2007.
//
// Usage:
//     REG = mser(IM[,DELTA[,MIN_SZ[,MAX_SZ]]]);
//
// Inputs:
//     IM        2D/3D image (real uint8)
//     DELTA     default 5, see Kristensen, 2007 (real double in [1, 127])
//     MIN_SZ    default 60 (real double >= 0)
//     MAX_SZ    default 14400 (real double >= 0)
//
// Outputs:
//     REG       binary masks for each MSER (logical)
//     COUNT     number of MSERs detected (real double)
//
// Compile:
//     mex -output mser mser/mser.cpp mexmser.cpp
//
// In case of problems with MacOS SDKs, consider:
//     https://xolotl.readthedocs.io/en/master/troubleshooting

#include "mex.h" // MATLAB
#include "matrix.h" // MATLAB
#include <vector>
#include <list>
#include <iostream>
#include <cmath> // pow
#include "mser/mser.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    using pixel_t = unsigned char;
    enum in_arg {IM=0, DELTA, MIN_SZ, MAX_SZ, NUM_IN};
    enum out_arg {MSERS=0, COUNT, NUM_OUT};
    // Validate number of arguments.
    if (nlhs > NUM_OUT)
    {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (nrhs > NUM_IN || nrhs < 1)
    {
        mexErrMsgTxt("Inputs are not: IM[,DELTA[,MIN_SZ[,MAX_SZ]]].");
    }
    // Validate image data.
    const mwSize num_dim = mxGetNumberOfDimensions(prhs[IM]);
    if (mxIsScalar(prhs[IM]) || !mxIsUint8(prhs[IM]) || mxIsComplex(prhs[IM])
        || num_dim<2 || num_dim>3)
    {
        mexErrMsgTxt("IM isn't a real 2D/3D array of type uint8.");
    }
    const mwSize *size = mxGetDimensions(prhs[IM]);
    const mwSize numx = size[0];
    const mwSize numy = size[1];
    const mwSize numz = num_dim>2 ? size[2] : 1;
    pixel_t *buf = (pixel_t *)mxGetData(prhs[IM]);
    std::vector<pixel_t> im(buf, buf + numx*numy*numz);
    // Validate options.
    int delta = 5, min_size = 60, max_size = 14400;
    if (nrhs > DELTA)
    {
        if (!mxIsScalar(prhs[DELTA]) || !mxIsDouble(prhs[DELTA])
            || mxIsComplex(prhs[DELTA]) || (double)mxGetScalar(prhs[DELTA]) < 1
            || (double)mxGetScalar(prhs[DELTA]) > pow(2, 7*sizeof(pixel_t))-1)
        {
            mexErrMsgTxt("DELTA isn't a double in [1, 127].");
        }
        delta = (int)mxGetScalar(prhs[DELTA]);
    }
    if (nrhs > MIN_SZ)
    {
        if (!mxIsScalar(prhs[MIN_SZ]) || !mxIsDouble(prhs[MIN_SZ])
            || mxIsComplex(prhs[MIN_SZ])
            || (double)mxGetScalar(prhs[MIN_SZ])<1)
        {
            mexErrMsgTxt("MIN_SZ isn't a double > 1.");
        }
        min_size = (int)mxGetScalar(prhs[MIN_SZ]);
    }
    if (nrhs > MAX_SZ)
    {
        if (!mxIsScalar(prhs[MAX_SZ]) || !mxIsDouble(prhs[MAX_SZ])
            || mxIsComplex(prhs[MAX_SZ]) || (double)mxGetScalar(prhs[MAX_SZ])<0
            || (double)mxGetScalar(prhs[MAX_SZ])<min_size)
        {
            mexErrMsgTxt("MAX_SZ isn't a double >= MIN_SZ.");
        }
        max_size = (int)mxGetScalar(prhs[MAX_SZ]);
    }
    // Detect MSERs.
    auto bins = bin_voxel_ind<pixel_t>(im);
    auto msers = detect_msers(bins, numx, numy, numz, delta, min_size,
        max_size);
    // Convert index lists to binary mask.
    const mwSize num_mser = msers.size();
    plhs[COUNT] = mxCreateDoubleScalar(num_mser);
    plhs[MSERS] = mxCreateCellMatrix(1, num_mser);
    int n = 0;
    for(const auto& ind_list : msers)
    {
        mxArray *dat = mxCreateLogicalArray(num_dim, size);
        mxLogical *mask = mxGetLogicals(dat);
        for (mwIndex j = 0; j < ind_list.size(); ++j)
        {
            mask[ ind_list[j] ] = true;
        }
        mxSetCell(plhs[MSERS], n++, dat);
    }
}
