/**
 * @file VectorizationPatternOption.h
 * @author L. Gall
 * @date 07/23/24
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
    inline namespace options {

    class VectorizationPatternOption : public Option<VectorizationPatternOption> {
        public:

            enum Value {

                p1xVec,

                p2xVecDiv2,

                pVecDiv2x2,

                pVecx1,

                pVecxVec
            };

            VectorizationPatternOption() = default;

            constexpr VectorizationPatternOption(Value option) : _value(option) {}

            constexpr operator Value() const { return _value; }

            static std::set<VectorizationPatternOption> getDiscouragedOptions() {
                return {};
            }

            static std::map<VectorizationPatternOption, std::string> getOptionNames() {
                return {
                    { VectorizationPatternOption::p1xVec, "1 x Vector Length"},
                    { VectorizationPatternOption::p2xVecDiv2, "2 x Vector Length / 2"},
                    { VectorizationPatternOption::pVecDiv2x2, "Vector Length / 2 x 2"},
                    { VectorizationPatternOption::pVecx1, "Vector Length x 1"},
                    { VectorizationPatternOption::pVecxVec, "Vector Length x Vector Length"},
                };
            }
        private:
            Value _value{Value(-1)};
    };

}
}