/*
 * Extract maximally stable extremal regions (MSERs) from a 2D or 3D image
 * (see Matas, J et al. Robust Wide Baseline Stereo from Maximally Stable
 * Extremal Regions. 2002). Implementation inspired by Kristensen, F et al.
 * Real-Time Extraction of Maximally Stable Extremal Regions on an FPGA.
 * 2007.
 *
 * Usage:
 *      REG = mser(IM[,DELTA[,MINSZ[,MAXSZ]]]);
 *
 * Inputs:
 *      IM          2D/3D image (real uint8)
 *      DELTA    	default 5, see Kristensen, 2007 (real double >= 1)
 *      MINSZ       default 60 (real double >= 0)
 *      MAXSZ       default 14400 (real double >= 0)
 *
 * Outputs:
 *      MSERS       binary masks for each MSER (logical)
 *      COUNT       number of MSERs detected (real double)
 *
 * Compile:
 *      mex mser.cpp
 */

#include "mex.h" // MATLAB
#include "matrix.h" // MATLAB
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <cassert>
#include <iostream>


// Efficient union find algorithm (see Kristensen, 2007).
class UnionFind {
private:
    size_t numRow;
    size_t numCol;
    size_t numSli;
    std::vector<int> regionMap;
    std::vector<size_t> nextPixel;
    std::unordered_set<size_t> refIndices;
    
private:
    void setU(size_t pixel, int value);
    void connect(size_t pixel);
    size_t findRef(size_t pixel) const;
    void swapLink(size_t ind1, size_t ind2);
    void mergeIf(size_t ref1, size_t ref2);
    size_t getNextPixel(size_t pixel) const;
    int getU(size_t pixel) const;

public:
    UnionFind(std::vector<size_t> dim);
    UnionFind(size_t numRow, size_t numCol=1, size_t numSli=1);
    void addPixel(size_t pixel);
    std::unordered_set<size_t>& getRefPoints();
    size_t getSize(size_t pixel) const;
    size_t getNumPixel() const;
    std::vector<size_t> getUnion(size_t pixel, size_t offset=0);
};


UnionFind::UnionFind(size_t numRow, size_t numCol, size_t numSli)
    :numRow( numRow )
    ,numCol( numCol )
    ,numSli( numSli )
    ,regionMap( getNumPixel(), 0 )
    ,nextPixel( getNumPixel(), 0 )
{
    assert(numRow > 1);
    assert((numCol>1 && numSli>1) || numSli==1);
    assert(getNumPixel() > 1);
}


void UnionFind::setU(size_t pixel, int value)
{
    regionMap[pixel] = value;
}


void UnionFind::connect(size_t pixel)
{
    nextPixel[pixel] = pixel;
    regionMap[pixel] = -1;
}


size_t UnionFind::findRef(size_t pixel) const
{
    assert(pixel < numRow*numCol*numSli);
    while (getU(pixel)>0)
    {
        pixel = getU(pixel);
    }
    return pixel;
}


void UnionFind::swapLink(size_t ind1, size_t ind2)
{
    size_t tmp = nextPixel[ind1];
    nextPixel[ind1] = nextPixel[ind2];
    nextPixel[ind2] = tmp;
}


void UnionFind::mergeIf(size_t ref1, size_t ref2)
{
    if (regionMap[ref1]==0 || regionMap[ref2]==0) // Ensure both placed!
    {
        return;
    }
    ref1 = findRef(ref1); // Ensure really have reference points.
    ref2 = findRef(ref2);
    if (ref1 == ref2)
    {
        return;
    }
    size_t size1 = getSize(ref1);
    size_t size2 = getSize(ref2);
    if(size2 < size1) // Insert smaller into larger cluster: stability.
    {
        size_t tmp = ref1;
        ref1 = ref2;
        ref2 = tmp;
        tmp = size1;
        size1 = size2;
        size2 = tmp;
    }
    setU(ref1, ref2); // Insert cluster #1 into #2.
    setU(ref2, -(size1+size2)); // Update new cluster size.
    swapLink(ref1, ref2);
    refIndices.emplace(ref2);
    if (size1 > 1) // Keep searching withing container to a minimum.
    {
        refIndices.erase(ref1); // Stop tracking smaller cluster.
    }
}


void UnionFind::addPixel(size_t pixel)
{
    assert(pixel < getNumPixel());
    connect(pixel);
    size_t row = pixel % numRow;
    size_t col = (pixel / numRow) % numCol;
    size_t sli = (pixel / numRow) / numCol;
    if (row > 0)
    {
        mergeIf(pixel, pixel-1);
    }
    if (col > 0)
    {
        mergeIf(pixel, pixel-numRow);
    }
    if (sli > 0)
    {
        mergeIf(pixel, pixel-numRow*numCol);
    }
    if (col < numCol-1)
    {
        mergeIf(pixel, pixel+numRow);
    }
    if (row < numRow-1)
    {
        mergeIf(pixel, pixel+1);
    }
    if (sli < numSli-1)
    {
        mergeIf(pixel, pixel+numRow*numCol);
    }
}


std::unordered_set<size_t>& UnionFind::getRefPoints()
{
    return refIndices;
}


size_t UnionFind::getSize(size_t pixel) const
{
    return -getU(findRef(pixel));
}


int UnionFind::getU(size_t pixel) const
{
    return regionMap[pixel];
}


size_t UnionFind::getNextPixel(size_t pixel) const
{
    return nextPixel[pixel];
}


size_t UnionFind::getNumPixel() const
{
    return numRow*numCol*numSli;
}


// Return current cluster, but skip the first offset pixels.
std::vector<size_t> UnionFind::getUnion(size_t pixel, size_t offset)
{
    size_t unionSize = getSize(pixel);
    size_t outSize = unionSize - offset;
    std::vector<size_t> outIndices(outSize);
    for (size_t i=0; i<unionSize; i++)
    {
        pixel = getNextPixel(pixel);
        if (i < offset)
        {
            continue;
        }
        outIndices[i-offset] = pixel;
    }
    return outIndices;
}


// Extract maximally stable extremal regions (MSERs) from a 2D or 3D image
// (see Matas, J et al. Robust Wide Baseline Stereo from Maximally Stable
// Extremal Regions. 2002). Implementation inspired by Kristensen, F et al.
// Real-Time Extraction of Maximally Stable Extremal Regions on an FPGA.
// 2007.
class MserDetector
{
public:
    using PixelType = unsigned char;
    static const size_t numIntLevel;
    static const size_t maxIntLevel;

private:
    UnionFind unionFind;
    std::vector< std::vector<size_t> > pixelsAtIntLevel;
    std::unordered_map<size_t, std::vector<size_t> > clusterSizes;
    std::unordered_map<size_t, double> growthRate;
    std::unordered_map<size_t, bool> lastRateLarger;

public:
    MserDetector(
        const std::vector<PixelType> im,
        size_t numRow,
        size_t numCol=1,
        size_t numSli=1);
    std::list< std::vector<size_t> > detect(
        size_t delta=5,
        size_t minSize=60,
        size_t maxSize=14400);
    std::vector<bool> createMask(std::vector<size_t> mser);

private:
    void reset();
    void addPixelsAtIntLevel(size_t intLevel);
    std::vector<size_t>& moveSizeWindow(size_t id, size_t size);
    bool checkMserAndUpdateRates(size_t id, double nextRate);
};

// Helper functions:

std::vector< std::vector<size_t> >
sortPixelInt(const std::vector<MserDetector::PixelType> im);

void invertImage(std::vector<MserDetector::PixelType>& im);

const size_t MserDetector::numIntLevel = pow(2, 8*sizeof(PixelType));
const size_t MserDetector::maxIntLevel = numIntLevel - 1;


MserDetector::MserDetector(
    const std::vector<PixelType> im,
    size_t numRow,
    size_t numCol,
    size_t numSli
)
    :unionFind( numRow, numCol, numSli )
    ,pixelsAtIntLevel( sortPixelInt(im) )
{
    assert(numRow > 1);
    assert((numCol>1 && numSli>1) || numSli==1);
    assert(im.size() > 1);
    assert(im.size() == numRow*numCol*numSli);
}


void MserDetector::addPixelsAtIntLevel(size_t intensity)
{
    std::vector<size_t>& pixelIndices = this->pixelsAtIntLevel[intensity];
    for (size_t index : pixelIndices)
    {
        this->unionFind.addPixel(index);
    }
}


std::vector<size_t>& MserDetector::moveSizeWindow(size_t id, size_t size)
{
    // Empty vector created if key not in unordered_map.
    std::vector<size_t>& window = this->clusterSizes[id];
    const bool isNew = window.size() == 0;
    if (isNew)
    {
        window.resize(size, 0);
        this->growthRate[id] = 0.0; // New entry into unordered_map.
        this->lastRateLarger[id] = false;
    }
    const size_t latestSize = unionFind.getSize(id);
    for (size_t i=0; i<size-1; i++)
    {
        window[i] = window[i+1];
    }
    window[size-1] = latestSize;
    return window;
}


void MserDetector::reset()
{
    this->clusterSizes.clear();
    this->growthRate.clear();
    this->lastRateLarger.clear();
}


bool MserDetector::checkMserAndUpdateRates(size_t id, double nextRate)
{
    bool isMser = lastRateLarger[id] && growthRate[id]<nextRate;
    this->lastRateLarger[id] = this->growthRate[id] > nextRate;
    this->growthRate[id] = nextRate;
    return isMser;
}


std::list< std::vector<size_t> > MserDetector::detect(
    size_t delta,
    size_t minSize,
    size_t maxSize)
{
    reset();
    assert(delta > 1);
    assert(minSize > 1);
    assert(maxSize > 1);
    std::list< std::vector<size_t> > output;
    const size_t windowSize = 2*delta+1;
    for (size_t level=0; level<MserDetector::numIntLevel; level++)
    {
        addPixelsAtIntLevel(level);
        std::unordered_set<size_t> clustIds = unionFind.getRefPoints();
        for (size_t id : clustIds)
        {
            std::vector<size_t> window = moveSizeWindow(id, windowSize);
            if (window[0] == 0) // Not full yet (size cannot decrease).
            {
                continue;
            }
            double nextRate = (window.back()-window[0])
                / (double)window[delta];
            bool isMser = checkMserAndUpdateRates(id, nextRate);
            size_t currSize = window[delta-1]; // Sizes one step ahead.
            if (!isMser || currSize<minSize || currSize>maxSize)
            {
                continue;
            }
            size_t offset = window.back()-currSize;
            std::vector<size_t> mser = unionFind.getUnion(id, offset);
            output.push_back(mser);
        }
    }
    return output;
}


std::vector<bool> MserDetector::createMask(std::vector<size_t> mser)
{
    size_t size = this->unionFind.getNumPixel();
    std::vector<bool> mask(size, false);
    std::vector<size_t>::iterator it;
    for(it = mser.begin(); it != mser.end(); ++it)
    {
        mask[*it] = true;
    }
    return mask;
}


std::vector<std::vector<size_t> > sortPixelInt(
    const std::vector<MserDetector::PixelType> im)
{
    // Number of intensity levels for image type.
    std::vector< std::vector<size_t> > ind(MserDetector::numIntLevel);
    for (size_t i=0; i<im.size(); i++)
    {
        MserDetector::PixelType pixelValue = im[i];
        ind[pixelValue].push_back(i);
    }
    return ind;
}


void invertImage(std::vector<MserDetector::PixelType>& im)
{
    for (auto& pixel : im)
    {
        pixel = MserDetector::maxIntLevel - pixel;
    }
}


/* MATLAB interface and input/output argument validation. */
void mexFunction(
        int nlhs,
        mxArray *plhs[],
        int nrhs,
        const mxArray *prhs[])
{
    enum inargs {IM=0, DELTA, MINSZ, MAXSZ, NUMIN};
    enum outargs {MSERS=0, COUNT, NUMOUT};
    
    /* Check number of input and output arguments. */
    if (nlhs > NUMOUT)
    {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (nrhs > NUMIN || nrhs < 1)
    {
        mexErrMsgTxt("Inputs aren't: IM[,DELTA[,MINSZ[,MAXSZ]]].");
    }
    
    /* Validate and access image. */
    const mwSize NDIM = mxGetNumberOfDimensions(prhs[IM]);
    if (mxIsScalar(prhs[IM]) || !mxIsUint8(prhs[IM]) ||
            mxIsComplex(prhs[IM]) || NDIM<2 || NDIM>3)
    {
        mexErrMsgTxt("IM isn't a 2D/3D array of type real uint8.");
    }
    const mwSize *size = mxGetDimensions(prhs[IM]);
    const mwSize numRow = size[0];
    const mwSize numCol = size[1];
    const mwSize numSli = NDIM>2 ? size[2] : 1;
    if (numRow<3 || numCol<3)
    {
        mexErrMsgTxt("IM is too small.");
    }
    using PixelType = MserDetector::PixelType;
    PixelType *im = (PixelType *)mxGetData(prhs[IM]);
    
    /* Validate options. */
    size_t delta=5, minSize=60, maxSize=14400;
    if (nrhs > DELTA)
    {
        if (!mxIsScalar(prhs[DELTA]) || !mxIsDouble(prhs[DELTA])
            || mxIsComplex(prhs[DELTA])
            || (double)mxGetScalar(prhs[DELTA])<1)
        {
            mexErrMsgTxt("DELTA isn't a scalar of type double >= 1.");
        }
        delta = (size_t)mxGetScalar(prhs[DELTA]);
    }
    if (nrhs > MINSZ)
    {
        if (!mxIsScalar(prhs[MINSZ]) || !mxIsDouble(prhs[MINSZ])
            || mxIsComplex(prhs[MINSZ])
            || (double)mxGetScalar(prhs[MINSZ])<0)
        {
            mexErrMsgTxt("MINSZ isn't a scalar of type double >=0.");
        }
        minSize = (size_t)mxGetScalar(prhs[MINSZ]);
    }
    if (nrhs > MAXSZ)
    {
        if (!mxIsScalar(prhs[MAXSZ]) || !mxIsDouble(prhs[MAXSZ])
            || mxIsComplex(prhs[MAXSZ])
            || (double)mxGetScalar(prhs[MAXSZ])<0)
        {
            mexErrMsgTxt("MAXSZ isn't a scalar of type double >=0.");
        }
        maxSize = (size_t)mxGetScalar(prhs[MAXSZ]);
    }
    if (minSize > maxSize)
    {
        mexErrMsgTxt("MINSZ is larger than MAXSZ.");
    }

    /* Detect MSERs. */
    const mwSize numVox = numRow * numCol * numSli;
    std::vector<PixelType> stdIm(im, im+numVox);
    MserDetector detector(stdIm, numRow, numCol, numSli);
    std::list< std::vector<size_t> > msers = detector.detect(delta,
        minSize, maxSize);

    /* Convert and prepare output arguments. */
    const mwSize numMser = msers.size();
    plhs[COUNT] = mxCreateDoubleScalar(numMser);
    plhs[MSERS] = mxCreateCellMatrix(1, numMser);
    size_t n = 0;
    for(std::vector<size_t> mserInd : msers)
    {
        const mwSize mserSize = mserInd.size();
        mxArray *dat = mxCreateLogicalArray(NDIM, size);
        mxLogical *mask = mxGetLogicals(dat);
        for (mwIndex j=0; j<mserSize; j++)
        {
            mask[ mserInd[j] ] = true;
        }
        mxSetCell(plhs[MSERS], n++, dat);
    }
}
