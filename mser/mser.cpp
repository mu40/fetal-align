#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <cassert>
#include <limits>
#include <algorithm>
#include "mser.hpp"

union_find::union_find(int dim_x, int dim_y, int dim_z)
    :num_x( dim_x )
    ,num_y( dim_y )
    ,num_z( dim_z )
    ,region_map( dim_x * dim_y * dim_z )
    ,next_pixel( dim_x * dim_y * dim_z, -1 )
{
    assert(dim_x > 0 && dim_y > 0 && dim_z > 0);
    assert(std::numeric_limits<int>::max() >= dim_x * dim_y * dim_z);
}

int union_find::find_ref(int pixel) const
{
    assert(next_pixel[pixel] > -1); // Must be placed.
    while (region_map[pixel] > -1)
    {
        pixel = region_map[pixel];
    }
    return pixel;
}

void union_find::merge_unions(int ind1, int ind2)
{
    if (next_pixel[ind1] < 0 || next_pixel[ind2] < 0) // Both placed?
    {
        return;
    }
    auto ref1 = find_ref(ind1);
    auto ref2 = find_ref(ind2);
    if (ref1 == ref2)
    {
        return;
    }
    auto size1 = -region_map[ref1];
    auto size2 = -region_map[ref2];
    if(size2 < size1) // Insert smaller into larger cluster for stability.
    {
        std::swap(ref1, ref2);
        std::swap(size1, size2);
    }
    region_map[ref1] = ref2; // Insert cluster 1 into 2.
    region_map[ref2] = -(size1 + size2); // Update size of cluster 2.
    std::swap(next_pixel[ref1], next_pixel[ref2]); // Swap links.
    ref_pts.emplace(ref2);
    if (size1 > 1) // Only tracked clusters of size N>1.
    {
        ref_pts.erase(ref1);
    }
}

void union_find::add_pixel(int ind)
{
    assert(ind < num_x * num_y * num_z);
    assert(next_pixel[ind] < 0); // Must be new.
    next_pixel[ind] = ind;
    region_map[ind] = -1;
    const auto x = ind % num_x;
    const auto y = (ind / num_x) % num_y;
    const auto z = (ind / num_x) / num_y;
    if (x > 0)
    {
        merge_unions(ind, ind - 1);
    }
    if (y > 0)
    {
        merge_unions(ind, ind - num_x);
    }
    if (z > 0)
    {
        merge_unions(ind, ind - num_x * num_y);
    }
    if (y < num_y-1)
    {
        merge_unions(ind, ind + num_x);
    }
    if (x < num_x-1)
    {
        merge_unions(ind, ind + 1);
    }
    if (z < num_z-1)
    {
        merge_unions(ind, ind + num_x * num_y);
    }
}

// Return the current set of reference points.
const std::unordered_set<int>& union_find::get_union_ids() const
{
    return ref_pts;
}

// Return size of union given its reference point. The function returns zero if
// the pixel is not a reference point or a single pixel only.
int union_find::get_size(int ref) const
{
    assert(ref < num_x * num_y * num_z);
    return (region_map[ref] < -1) ? -region_map[ref] : 0;
}

// Return pixel indices of a union, skipping N=offset. An empty container is
// returned if the index is not a reference point or a single pixel.
std::vector<int> union_find::get_pixel_ind(int ref, int offset=0) const
{
    const int union_size = get_size(ref); // Zero if not ref or in range.
    const int out_size = union_size - offset;
    std::vector<int> out(out_size);
    for (int i=0; i<union_size; i++)
    {
        ref = next_pixel[ref]; // Ref point added last.
        if (i < offset)
        {
            continue;
        }
        out[i-offset] = ref;
    }
    return out;
}

// Extract maximally stable extremal regions (MSERs) from a 2D or 3D image
// (see Matas, J et al. Robust Wide Baseline Stereo from Maximally Stable
// Extremal Regions. 2002). Implementation inspired by Kristensen, F et al.
// Real-Time Extraction of Maximally Stable Extremal Regions on an FPGA. 2007.
std::list< std::vector<int> > detect_msers(
    const std::vector< std::vector<int> > binned_vox_ind,
    int num_x,
    int num_y,
    int num_z,
    int delta,
    int min_size,
    int max_size)
{
    assert(num_x > 0 && num_y > 0 && num_z > 0);
    assert(delta > 0 && delta < static_cast<int>(binned_vox_ind.size()/2));
    assert(min_size > 1);
    assert(max_size >= min_size);
    std::unordered_map< int, std::vector<int> > cluster_sizes;
    std::unordered_map< int, double > rate, last;
    union_find clusters(num_x, num_y, num_z);
    std::list< std::vector<int> > out;
    const auto win_size = 2 * delta + 1;
    for (const auto& bin : binned_vox_ind) // Iterate intensity levels.
    {
        for (auto ind : bin)
        {
            clusters.add_pixel(ind);
        }
        for (auto id : clusters.get_union_ids())
        {
            // Update sliding cluster-size window; don't loop if empty.
            std::vector<int>& win = cluster_sizes[id];
            for (int i = 0; i < static_cast<int>(win.size()) - 1; ++i)
            {
                win[i] = win[i+1];
            }
            if (win.size() == 0)
            {
                win.resize(win_size, 0);
                rate[id] = last[id] = 0.0; // Unknown key creates entry.
            }
            win.back() = clusters.get_size(id);
            // Check if MSER and update growth rates once is full enough.
            if (win[delta] == 0)
            {
                continue;
            }
            const auto next = (win.back() - win[0]) /
                static_cast<double>(win[delta]);
            const bool is_min = last[id] > rate[id] && rate[id] < next;
            last[id] = rate[id];
            rate[id] = next;
            const auto mser_size = win[delta-1]; // Sizes one step ahead.
            if (!is_min || mser_size < min_size || mser_size > max_size)
            {
                continue;
            }
            const auto offset = win.back() - mser_size;
            const std::vector<int> mser = clusters.get_pixel_ind(id, offset);
            out.push_back(mser);
        }
    }
    return out;
}
