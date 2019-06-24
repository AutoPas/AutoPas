// Copyright 2016, Tobias Hermann.
// https://github.com/Dobiasd/frugally-deep
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "fdeep/layers/pooling_2d_layer.hpp"

#include <algorithm>
#include <limits>
#include <string>

namespace fdeep { namespace internal
{

FDEEP_FORCE_INLINE tensor5 max_pool_2d(
    std::size_t pool_height, std::size_t pool_width,
    std::size_t strides_y, std::size_t strides_x,
    bool channels_first,
    padding pad_type,
    bool use_offset,
    const tensor5& in)
{
    const float_type invalid = std::numeric_limits<float_type>::lowest();

    const std::size_t feature_count = channels_first
        ? in.shape().height_
        : in.shape().depth_
        ;

    const std::size_t in_height = channels_first
        ? in.shape().width_
        : in.shape().height_
        ;

    const std::size_t in_width = channels_first
        ? in.shape().depth_
        : in.shape().width_
        ;

    const auto conv_cfg = preprocess_convolution(
        shape2(pool_height, pool_width),
        shape2(strides_y, strides_x),
        pad_type, use_offset, in_height, in_width);

    int pad_top_int = static_cast<int>(conv_cfg.pad_top_);
    int pad_left_int = static_cast<int>(conv_cfg.pad_left_);
    const std::size_t offset_y = conv_cfg.offset_y_;
    const std::size_t offset_x = conv_cfg.offset_x_;
    const std::size_t out_height = conv_cfg.out_height_;
    const std::size_t out_width = conv_cfg.out_width_;

    if (channels_first)
    {
        tensor5 out(shape5(1, 1, feature_count, out_height, out_width), 0);

        for (std::size_t z = 0; z < feature_count; ++z)
        {
            for (std::size_t y = 0; y < out_height; ++y)
            {
                for (std::size_t x = 0; x < out_width; ++x)
                {
                    float_type val = std::numeric_limits<float_type>::lowest();
                    for (std::size_t yf = 0; yf < pool_height; ++yf)
                    {
                        int in_get_y = static_cast<int>(offset_y + strides_y * y + yf) - pad_top_int;
                        for (std::size_t xf = 0; xf < pool_width; ++xf)
                        {
                            int in_get_x = static_cast<int>(offset_x + strides_x * x + xf) - pad_left_int;
                            const auto current = in.get_x_z_padded(invalid, z, in_get_y, in_get_x);
                            val = std::max(val, current);
                        }
                    }

                    out.set(0, 0, z, y, x, val);
                }
            }
        }
        return out;
    }
    else
    {
        tensor5 out(shape5(1, 1, out_height, out_width, feature_count), 0);

        for (std::size_t y = 0; y < out_height; ++y)
        {
            for (std::size_t x = 0; x < out_width; ++x)
            {
                for (std::size_t z = 0; z < feature_count; ++z)
                {
                    float_type val = std::numeric_limits<float_type>::lowest();
                    for (std::size_t yf = 0; yf < pool_height; ++yf)
                    {
                        int in_get_y = static_cast<int>(offset_y + strides_y * y + yf) - pad_top_int;
                        for (std::size_t xf = 0; xf < pool_width; ++xf)
                        {
                            int in_get_x = static_cast<int>(offset_x + strides_x * x + xf) - pad_left_int;
                            const auto current = in.get_y_x_padded(invalid, in_get_y, in_get_x, z);
                            val = std::max(val, current);
                        }
                    }

                    out.set(0, 0, y, x, z, val);
                }
            }
        }
        return out;
    }
}

class max_pooling_2d_layer : public pooling_2d_layer
{
public:
    explicit max_pooling_2d_layer(const std::string& name,
        const shape2& pool_size, const shape2& strides, bool channels_first, padding p,
        bool padding_valid_uses_offset, bool padding_same_uses_offset) :
        pooling_2d_layer(name, pool_size, strides, channels_first, p,
            padding_valid_uses_offset, padding_same_uses_offset)
    {
    }
protected:
    tensor5 pool(const tensor5& in) const override
    {
        if (pool_size_ == shape2(2, 2) && strides_ == shape2(2, 2))
            return max_pool_2d(2, 2, 2, 2, channels_first_, padding_, use_offset(), in);
        else if (pool_size_ == shape2(4, 4) && strides_ == shape2(4, 4))
            return max_pool_2d(4, 4, 4, 4, channels_first_, padding_, use_offset(), in);
        else
            return max_pool_2d(
                pool_size_.height_, pool_size_.width_,
                strides_.height_, strides_.width_,
                channels_first_, padding_, use_offset(), in);
    }
};

} } // namespace fdeep, namespace internal
