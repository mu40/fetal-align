#ifndef MSER_HPP
#define MSER_HPP

#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <cassert>
#include <type_traits>
#include <limits>

// Efficient union-find algorithm inspired by Kristensen (2007) using a region
// map U of the same size and dimension as the image. Pixels with U>=0 are part
// of the same union (or cluster) as pixel U. Pixels with U<0 are the single
// reference point of a union consisting of -U pixels. To efficiently find the
// pixels of a union, a circular list of links holds the position L>=0 of the
// next pixel in the union. A value of L<0 means that a pixel has not been
// placed in the region map yet. We also keep track of all unions (here defined
// as consisting of N>1 pixels) for efficient access.
class union_find
{
private:
    int num_x;
    int num_y;
    int num_z;
    std::vector<int> region_map; // U.
    std::vector<int> next_pixel; // L.
    std::unordered_set<int> ref_pts;

private:
    int find_ref(int pixel) const;
    void merge_unions(int ref1, int ref2);

public:
    union_find(int num_x, int num_y=1, int num_z=1);
    void add_pixel(int ind);
    const std::unordered_set<int>& get_union_ids() const;
    int get_size(int ref) const;
    std::vector<int> get_pixel_ind(int ref, int offset) const;
};


std::list< std::vector<int> > detect_msers(
    const std::vector< std::vector<int> > binned_vox_ind,
    int num_x,
    int num_y,
    int num_z,
    int delta=5,
    int min_size=60,
    int max_size=14400);

// Create binary mask of pixel type T from list of indices.
template<class T>
std::vector<T> fill_mask(
    const std::vector<int>& ind,
    int num_mask_vox,
    T bg=0,
    T fg=1)
{
    std::vector<T> mask = std::vector<T>(num_mask_vox, bg);
    for(auto i : ind)
    {
        mask.at(i) = fg;
    }
    return mask;
}

// Bin voxel indices of an image of integral type according to intensity values.
template<class T>
std::vector< std::vector<int> > bin_voxel_ind(const std::vector<T>& im)
{
    assert(std::is_integral<T>::value);
    const T min = std::numeric_limits<T>::lowest();
    const T max = std::numeric_limits<T>::max();
    std::vector< std::vector<int> > bins(max - min + 1);
    for (int i = 0; i < im.size(); i++)
    {
        const int val = im[i] - min;
        bins[val].push_back(i);
    }
    return bins;
}

#ifdef USE_ITK
#include <itkImage.h>
#include "itk.hpp"

// Bin voxel indices of an image of integral type according to intensity values.
template<class TImage>
std::vector< std::vector<int> >
bin_voxel_ind(const typename TImage::Pointer im)
{
    using pixel_t = typename TImage::PixelType;
    using iter_t = itk::ImageRegionConstIterator<TImage>;
    assert(std::is_integral<pixel_t>::value);
    assert(im);
    iter_t it(im, im->GetLargestPossibleRegion());
    const pixel_t min = std::numeric_limits<pixel_t>::lowest();
    const pixel_t max = std::numeric_limits<pixel_t>::max();
    std::vector< std::vector<int> > bins(max - min + 1);
    int i = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        const int val = it.Get() - min;
        bins[val].push_back(i++);
    }
    return bins;
}

// Create binary mask from list of indices. The function overwrites an existing
// image buffer or allocates a new image using information from a reference.
template<class TImage, class TImageRef=TImage>
typename TImage::Pointer fill_mask(
    const std::vector<int>& ind,
    typename TImage::Pointer mask,
    const typename TImageRef::Pointer meta_info=nullptr,
    typename TImage::PixelType bg=0,
    typename TImage::PixelType fg=1)
{
    assert(meta_info || mask);
    if (mask)
    {
        mask->Update();
    }
    else
    {
        meta_info->Update();
        mask = clone_image<TImageRef, TImage>(meta_info);
    }
    mask->FillBuffer( bg );
    typename TImage::PixelType* buf = mask->GetBufferPointer();
    const auto num_vox = count_voxels<TImage>(mask);
    for (auto i : ind)
    {
        if (i < num_vox) {
            buf[i] = fg;
        }
    }
    return mask;
}

#endif // USE_ITK

#endif // MSER_HPP
