// Copyright 2016, Tobias Hermann.
// https://github.com/Dobiasd/frugally-deep
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "fdeep/common.hpp"

#include <cstddef>
#include <cstdlib>
#include <string>

namespace fdeep { namespace internal
{

class tensor5_pos
{
public:
    // The dimensions are right-aligned (left-padded) compared to Keras.
    // I.e., if you have a position (or shape) of (a, b) in Keras
    // it corresponds to (0, 0, 0, a, b) in frugally-deep.
    explicit tensor5_pos(
        std::size_t pos_dim_5,
        std::size_t pos_dim_4,
        std::size_t y,
        std::size_t x,
        std::size_t z) :
            pos_dim_5_(pos_dim_5),
            pos_dim_4_(pos_dim_4),
            y_(y),
            x_(x),
            z_(z)
    {
    }

    std::size_t pos_dim_5_;
    std::size_t pos_dim_4_;
    std::size_t y_;
    std::size_t x_;
    std::size_t z_;
};

inline std::size_t get_tensor5_pos_dimension_by_index(const tensor5_pos& s,
    const std::size_t idx)
{
    if (idx == 0)
        return s.pos_dim_5_;
    if (idx == 1)
        return s.pos_dim_4_;
    if (idx == 2)
        return s.y_;
    if (idx == 3)
        return s.x_;
    if (idx == 4)
        return s.z_;
    raise_error("Invalid tensor5_pos index.");
    return 0;
}

inline tensor5_pos change_tensor5_pos_dimension_by_index(const tensor5_pos& in,
    const std::size_t idx, const std::size_t dim)
{
    tensor5_pos out = in;
    if (idx == 0)
        out.pos_dim_5_ = dim;
    else if (idx == 1)
        out.pos_dim_4_ = dim;
    else if (idx == 2)
        out.y_ = dim;
    else if (idx == 3)
        out.x_ = dim;
    else if (idx == 4)
        out.z_ = dim;
    else
        raise_error("Invalid tensor5_pos index.");
    return out;
}

} } // namespace fdeep, namespace internal
